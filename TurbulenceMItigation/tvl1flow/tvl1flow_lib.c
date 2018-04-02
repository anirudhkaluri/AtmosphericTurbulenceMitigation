


#ifndef DUAL_TVL1_OPTIC_FLOW_H
#define DUAL_TVL1_OPTIC_FLOW_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "mask.c"
#include "warp.c"
#include "zoom.c"

#define MAX_ITERATIONS 150




static void *xmalloc(size_t size)
{
    void *r = malloc(size);
    if (!r)
    {
        double sm = size / (0x100000 * 1.0);
        fprintf(stderr, "malloc: out of memory when requeting "
                "%zu bytes (%gMB)\n", size, sm);
        abort();
    }
    return r;
}




void Dual_TVL1_optic_flow(
    float *I0,                   //source image
    float *I1,                   //target image
    float *u1,                   //x component of the optical flow
    float *u2,                   //y component of the optical flow
    const int   nx,              //image width
    const int   ny,              //image height
    const float tau    ,   //time step
    const float lambda ,   //weight parameter for the data term
    const float theta   ,    //weight parameter for (u - v)²
    const int   warps  ,      //number of warpings per scale
    const float epsilon    , //tolerance for numerical convergence
    const bool  verbose   //enable/disable the verbose mode
)
{
    const int   size = nx * ny;
    const float l_t = lambda * theta;

    size_t sf = sizeof(float);
    float *I1x    = (float*)xmalloc(size*sf);
    float *I1y    = (float*)xmalloc(size*sf);
    float *I1w    = (float*)xmalloc(size*sf);
    float *I1wx   = (float*)xmalloc(size*sf);
    float *I1wy   = (float*)xmalloc(size*sf);
    float *rho_c  = (float*)xmalloc(size*sf);
    float *v1     = (float*)xmalloc(size*sf);
    float *v2     = (float*)xmalloc(size*sf);
    float *p11    = (float*)xmalloc(size*sf);
    float *p12    = (float*)xmalloc(size*sf);
    float *p21    = (float*)xmalloc(size*sf);
    float *p22    = (float*)xmalloc(size*sf);
    float *div    = (float*)xmalloc(size*sf);
    float *grad   = (float*)xmalloc(size*sf);
    float *div_p1 = (float*)xmalloc(size*sf);
    float *div_p2 = (float*)xmalloc(size*sf);
    float *u1x    = (float*)xmalloc(size*sf);
    float *u1y    = (float*)xmalloc(size*sf);
    float *u2x    = (float*)xmalloc(size*sf);
    float *u2y    = (float*)xmalloc(size*sf);
    float *mask   = (float*)xmalloc(size*sf);

    centered_gradient(I1, I1x, I1y, nx, ny);


    //initialization of p
    for(int i = 0; i < size; i++)
    {
        p11[i] = p12[i] = 0.0;
        p21[i] = p22[i] = 0.0;
    }


    for(int warpings = 0; warpings < warps; warpings++)
    {
        //compute the warping of the target image and its derivatives
        bilinear_interpolation(I1,  u1, u2, I1w,  nx, ny, mask);
        bilinear_interpolation(I1x, u1, u2, I1wx, nx, ny, mask);
        bilinear_interpolation(I1y, u1, u2, I1wy, nx, ny, mask);

        for(int i = 0; i < size; i++)
        {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            //store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            //compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }


        int n = 0;
        float error = 1000;
        while(error > epsilon && n < MAX_ITERATIONS)
        {
            n++;
            //estimate the values of the variable (v1, v2)
            for(int i = 0; i < size; i++)
            {
                const float rho = rho_c[i] + (I1wx[i] * u1[i]
                                                  + I1wy[i] * u2[i]);

                float d1, d2;

                if(rho < - l_t * grad[i])
                {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                }
                else
                {
                    if(rho > l_t * grad[i])
                    {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    }
                    else
                    {
                        if(grad[i] < 1E-10)
                            d1 = d2 = 0;
                        else
                        {
                            d1 = - rho * I1wx[i] / grad[i];
                            d2 = - rho * I1wy[i] / grad[i];
                        }
                    }
                }

                v1[i] = u1[i] + mask[i] * d1;
                v2[i] = u2[i] + mask[i] * d2;
            }


            //compute the divergence of the dual variable (p1, p2)
            divergence(p11, p12, div_p1, nx ,ny);
            divergence(p21, p22, div_p2, nx ,ny);

            //estimate the values of the optical flow (u1, u2)
            error = 0.0;
            for(int i = 0; i < size; i++)
            {
                const float u1k = u1[i];
                const float u2k = u2[i];

                u1[i] = v1[i] + theta * div_p1[i];
                u2[i] = v2[i] + theta * div_p2[i];

                error += (u1[i] - u1k) * (u1[i] - u1k) +
                            (u2[i] - u2k) * (u2[i] - u2k);
            }

            error /= size;

            //compute the gradient of the optical flow (Du1, Du2)
            forward_gradient(u1, u1x, u1y, nx ,ny);
            forward_gradient(u2, u2x, u2y, nx ,ny);

            //estimate the values of the dual variable (p1, p2)
            for(int i = 0; i < size; i++)
            {
                const float taut = tau / theta;
                const float g1   = sqrt(u1x[i] * u1x[i] + u1y[i] * u1y[i]);
                const float g2   = sqrt(u2x[i] * u2x[i] + u2y[i] * u2y[i]);
                const float ng1  = 1.0 + taut * g1;
                const float ng2  = 1.0 + taut * g2;

                p11[i] = (p11[i] + taut * u1x[i]) / ng1;
                p12[i] = (p12[i] + taut * u1y[i]) / ng1;
                p21[i] = (p21[i] + taut * u2x[i]) / ng2;
                p22[i] = (p22[i] + taut * u2y[i]) / ng2;
            }
        }
    }

    free(I1x);
    free(I1y);
    free(I1w);
    free(I1wx);
    free(I1wy);
    free(rho_c);
    free(v1);
    free(v2);
    free(p11);
    free(p12);
    free(p21);
    free(p22);
    free(div);
    free(grad);
    free(div_p1);
    free(div_p2);
    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);
    free(mask);
}

static void getminmax(float *min, float *max, const float *x, int n)
{
    *min = *max = x[0];
    for (int i = 1; i < n; i++)
    {
        if (x[i] < *min)
            *min = x[i];
        if (x[i] > *max)
            *max = x[i];
    }
}

/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization(
    const float *I0,  //input image0
    const float *I1,  //input image1
    float *I0n,       //normalized output image0
    float *I1n,       //normalized output image1
    int size          //size of the image
)
{
    float max0, max1, min0, min1;
    getminmax(&min0, &max0, I0, size);
    getminmax(&min1, &max1, I1, size);

    const float max = max0 > max1 ? max0 : max1;//std::max(max0, max1);
    const float min = min0 < min1 ? min0 : min1;//std::min(min0, min1);
    const float den = max - min;

    if(den > 0)

        for(int i = 0; i < size; i++)
        {
            I0n[i] = 255.0 * (I0[i] - min) / den;
            I1n[i] = 255.0 * (I1[i] - min) / den;
        }

    else

        for(int i = 0; i < size; i++)
        {
            I0n[i] = I0[i];
            I1n[i] = I1[i];
        }
}


/**
 *
 * Function to compute the optical flow using multiple scales
 *
 **/
void Dual_TVL1_optic_flow_multiscale(
    float *I0,                   //source image
    float *I1,                   //target image
    float *u1,                   //x component of the optical flow
    float *u2,                   //y component of the optical flow
    const int   nxx,             //image width
    const int   nyy,             //image height
    const float tau     ,  //time step
    const float lambda  ,  //weight parameter for the data term
    const float theta    ,   //weight parameter for (u - v)²
    const int   nscales ,     //number of scales
    const float zfactor ,   //reduce factor for creating the image piramid
    const int   warps   ,     //number of warpings per scale
    const float epsilon    , //tolerance for numerical convergence
    const bool  verbose   //enable/disable the verbose mode
)
{
    int size = nxx * nyy;
    float **I0s = (float**)xmalloc(nscales * sizeof(float*));
    float **I1s = (float**)xmalloc(nscales * sizeof(float*));
    float **u1s = (float**)xmalloc(nscales * sizeof(float*));
    float **u2s = (float**)xmalloc(nscales * sizeof(float*));
    int *nx = (int*)xmalloc(nscales * sizeof(int));
    int *ny = (int*)xmalloc(nscales * sizeof(int));

    I0s[0] = (float*)xmalloc(size*sizeof(float));
    I1s[0] = (float*)xmalloc(size*sizeof(float));

    u1s[0] = u1;
    u2s[0] = u2;
    nx [0] = nxx;
    ny [0] = nyy;

    //normalize the images between 0 and 255
    image_normalization(I0, I1, I0s[0], I1s[0], size);

    //initialize the flow at the finest scale
    for(int i = 0; i < size; i++)

        u1s[0][i] = u2s[0][i] = 0.0;

    //create the scales
    for(int s = 1; s < nscales; s++)
    {
        zoom_size(nx[s-1], ny[s-1], &nx[s], &ny[s], zfactor);
        const int sizes = nx[s] * ny[s];

        I0s[s] = (float*)xmalloc(sizes*sizeof(float));
        I1s[s] = (float*)xmalloc(sizes*sizeof(float));
        u1s[s] = (float*)xmalloc(sizes*sizeof(float));
        u2s[s] = (float*)xmalloc(sizes*sizeof(float));

        if(s == nscales - 1)
            //initialize the flow at the coarsest scale
            for(int i = 0; i < sizes; i++)
                u1s[s][i] = u2s[s][i] = 0.0;

        //zoom in the images to create the pyramidal structure
        zoom(I0s[s-1], I0s[s], nx[s-1], ny[s-1], zfactor);
        zoom(I1s[s-1], I1s[s], nx[s-1], ny[s-1], zfactor);
    }

    for(int s = nscales-1; s >= 0; s--)
    {
        //compute the optical flow at the current scale
        Dual_TVL1_optic_flow(
            I0s[s], I1s[s], u1s[s], u2s[s], nx[s], ny[s],
            tau, lambda, theta, warps, epsilon, verbose
        );

        if(s)
        {
            zoom_out(u1s[s], u1s[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
            zoom_out(u2s[s], u2s[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);

            int size = nx[s-1] * ny[s-1];

            for(int i = 0; i < size; i++)
            {
                u1s[s-1][i] *= (float) 1.0 / zfactor;
                u2s[s-1][i] *= (float) 1.0 / zfactor;
            }
        }
    }

    for(int i = 1; i < nscales; i++)
    {
        free(I0s[i]);
        free(I1s[i]);
        free(u1s[i]);
        free(u2s[i]);
    }
    free(I0s[0]);
    free(I1s[0]);

    free(I0s);
    free(I1s);
    free(u1s);
    free(u2s);
    free(nx);
    free(ny);
}


#endif
