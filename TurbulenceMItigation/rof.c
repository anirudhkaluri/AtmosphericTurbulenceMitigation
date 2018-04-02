

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "iio.h"

static void *rof_xmalloc(size_t size)
{
    void *r = malloc(size);
    if (!r)
    {
        double sm = size / (0x100000 * 1.0);
        fprintf(stderr, "rof_xmalloc: out of memory when requeting "
                "%zu bytes (%gMB)\n", size, sm);
        abort();
    }
    return r;
}


static void rof_divergence(
    const float *v1, //x component of the vector field
    const float *v2, //y component of the vector field
    float *div,      //output divergence
    const int nx,    //image width
    const int ny     //image height
)
{
				//apply the divergence to the center body of the image
        for(int i = 1; i < ny-1; i++)
        {
            for(int j = 1; j < nx-1; j++)
            {
                const int p  = i * nx + j;
                const int p1 = p - 1;
                const int p2 = p - nx;

                const float v1x = v1[p] - v1[p1];
                const float v2y = v2[p] - v2[p2];

                div[p] = v1x + v2y;
            }
        }

        for(int j = 1; j < nx-1; j++)
        {
            const int p = (ny-1) * nx + j;

            div[j] = v1[j] - v1[j-1] + v2[j];
            div[p] = v1[p] - v1[p-1] - v2[p-nx];
        }

        //apply the divergence to the first and last columns
        for(int i = 1; i < ny-1; i++)
        {
            const int p1 = i * nx;
            const int p2 = (i+1) * nx - 1;

            div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
            div[p2] = -v1[p2-1] + v2[p2] - v2[p2 - nx];

        }

        div[0]         =  v1[0] + v2[0];
        div[nx-1]      = -v1[nx - 2] + v2[nx - 1];
        div[(ny-1)*nx] =  v1[(ny-1)*nx] - v2[(ny-2)*nx];
        div[ny*nx-1]   = -v1[ny*nx - 2] - v2[(ny-1)*nx - 1];
}


/**
 *
 * Function to compute the gradient with forward differences
 * Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
 *
 **/
static void rof_forward_gradient(
    const float *f, //input image
    float *fx,      //computed x derivative
    float *fy,      //computed y derivative
    float *grad_norm, // gradient norm output
    const int nx,   //image width
    const int ny    //image height
)
{
        //apply the gradient to the center body of the image
        for(int i = 0; i < ny-1; i++)
        {
            for(int j = 0; j < nx-1; j++)
            {
                const int p  = i * nx + j;
                const int p1 = p + 1;
                const int p2 = p + nx;

                fx[p] = f[p1] - f[p];
                fy[p] = f[p2] - f[p];
								grad_norm[p] = sqrt(powf(fx[p],2)+powf(fy[p],2));
            }
        }

        //apply the gradient to the last row
        for(int j = 0; j < nx-1; j++)
        {
            const int p = (ny-1) * nx + j;

            fx[p] = f[p+1] - f[p];
            fy[p] = 0;
						grad_norm[p] = sqrt(powf(fx[p],2)+powf(fy[p],2));
        }

        //apply the gradient to the last column
        for(int i = 1; i < ny; i++)
        {
            const int p = i * nx-1;

            fx[p] = 0;
            fy[p] = f[p+nx] - f[p];
						grad_norm[p] = sqrt(powf(fx[p],2)+powf(fy[p],2));
        }

        fx[ny * nx - 1] = 0;
        fy[ny * nx - 1] = 0;
				grad_norm[ny * nx - 1] = sqrt(powf(fx[ny * nx - 1],2)+powf(fy[ny * nx - 1],2));
}



void RudinOsherFatemi(
    float *in, // input image
    float *d, // denoised output
    int w, // image width
    int h, // image height
    float lambda, // regularization parameter
    int itemax // number max of iterations
)
{

    float *dxw = NULL;
    float *dyw = NULL;
    float *norm = NULL;
    float *phix = NULL;
    float *phiy = NULL;

    int i,j,b;
    int it = 0;
    float tau = 0.12;
    float erreur = 0.001; //to test the convergence
    float err,maxerr;
    float npx,npy,co;

    dxw = rof_xmalloc(sizeof(float)*w*h);
    dyw = rof_xmalloc(sizeof(float)*w*h);
    norm = rof_xmalloc(sizeof(float)*w*h);
    phix = rof_xmalloc(sizeof(float)*w*h);
    phiy = rof_xmalloc(sizeof(float)*w*h);

    //initialization of d,phix,phiy to 0
    for (i=0; i<w; i++)
    {
        for (j=0; j<h; j++)
        {
            d[i+j*w] = 0.0;
            phix[i+j*w] = 0.0;
            phiy[i+j*w] = 0.0;
        }
    }

    co = tau/lambda;

    do
    {
        //we compute d=lambda*d-in
        for (i=0; i<w; i++)
        {
            for (j=0; j<h; j++)
            {
                d[i+j*w] = lambda * d[i+j*w] - in[i+j*w];
            }
        }

        //we compute the derivative and gradient norm of d
				rof_forward_gradient(in,dxw,dyw,norm,w,h);

        //we update the projector
        maxerr = 0;

        for (i=0; i<w; i++)
        {
            for (j=0; j<h; j++)
            {
                npx = (phix[i+j*w] + co*dxw[i+j*w])/(1. + co*norm[i+j*w]);
                npy = (phiy[i+j*w] + co*dyw[i+j*w])/(1. + co*norm[i+j*w]);
                err = sqrt(powf(phix[i+j*w]-npx,2)
                            + powf(phiy[i+j*w] - npy,2));
                if (err > maxerr) maxerr = err;
                phix[i+j*w] = npx;
                phiy[i+j*w] = npy;
            }
        }

        rof_divergence(phix,phiy,d,w,h);

        it++;
    }
    while ((maxerr > erreur) && (it < itemax));

    for (i=0; i<w; i++)
    {
        for (j=0; j<h; j++)
        {
            d[i+j*w] = in[i+j*w] - lambda * d[i+j*w];
        }
    }

    free(dxw);
    free(dyw);
    free(norm);
    free(phix);
    free(phiy);

    return;
}
