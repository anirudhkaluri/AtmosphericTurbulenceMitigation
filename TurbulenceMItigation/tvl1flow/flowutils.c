

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

#include "tvl1flow_lib.c"
#include "getpixel.c"
#include "bicubic.c"

#include "smapa.h"


static float evaluate_bilinear_cell(float a, float b, float c, float d,
                                    float x, float y)
{
    float r = 0;
    r += a * (1-x) * (1-y);
    r += b * ( x ) * (1-y);
    r += c * (1-x) * ( y );
    r += d * ( x ) * ( y );
    return r;
}

static float getsample(float *fx, int w, int h, int pd, int i, int j, int l)
{
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return 0;
    float (*x)[w][pd] = (void*)fx;
    return x[j][i][l];
}

static float getsamplen(float *fx, int w, int h, int pd, int i, int j, int l)
{
    if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
        return NAN;
    float (*x)[w][pd] = (void*)fx;
    return x[j][i][l];
}

static void bilinear_interpolation_at(float *result,
                                      float *x, int w, int h, int pd,
                                      float p, float q)
{
    int ip = p;
    int iq = q;
    FORL(pd)
    {
        float a = getsamplen(x, w, h, pd, ip  , iq  , l);
        float b = getsamplen(x, w, h, pd, ip+1, iq  , l);
        float c = getsamplen(x, w, h, pd, ip  , iq+1, l);
        float d = getsamplen(x, w, h, pd, ip+1, iq+1, l);
        float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
        result[l] = r;
    }
}

SMART_PARAMETER_SILENT(BACKDIV,0)

static void compute_flow_div(float *d, float *u, int w, int h)
{
    getsample_operator p = getsample_1;
    FORJ(h) FORI(w)
    {
        float ux = 0.5*(p(u,w,h,2, i+1, j, 0) - p(u,w,h,2, i-1, j, 0));
        float vy = 0.5*(p(u,w,h,2, i, j+1, 1) - p(u,w,h,2, i, j-1, 1));
        d[j*w+i] = ux + vy;
    }
}

/**
 *
 * Function to compute the warping of the input image by the
 * provided optical flow
 *
 **/
// Apply the deformation field flo to image in
void warpwithflow(float *ou, float *flo, float *pin, int w, int h, int pd)
{
    float (*out)[w][pd] = (void*)ou;
    float (*in)[w][pd] = (void*)pin;
    float (*flow)[w][2] = (void*)flo;
    float *flowdiv = NULL;

    if (BACKDIV() > 0)
    {
        flowdiv = xmalloc(w * h * sizeof*flowdiv);
        compute_flow_div(flowdiv, flo, w, h);
    }

    FORJ(h) FORI(w)
    {
        float p[2] = {i + flow[j][i][0], j + flow[j][i][1]};
        float result[pd];
        bicubic_interpolation(result, pin, w, h, pd, p[0], p[1]);
        float factor = 1;
        if (flowdiv)
            factor = exp(BACKDIV() * flowdiv[j*w+i]);
        FORL(pd)
        out[j][i][l] = factor * result[l];
    }

    if (flowdiv) free(flowdiv);
}

static bool checkbounds(int a, int x, int b)
{
    return a <= x && x < b;
}

static void flowinv_init(float *v, float *u, int w, int h)
{
    float (*U)[w][2] = (void*)u;
    float (*V)[w][2] = (void*)v;

    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
            for (int l = 0; l < 2; l++)
                V[j][i][l] = -U[j][i][l];

    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            float px = i + U[j][i][0];
            float py = j + U[j][i][1];
            int ipx = px;
            int ipy = py;
            if (checkbounds(0, ipx, w) && checkbounds(0, ipy, h))
                for (int l = 0; l < 2; l++)
                    V[ipy][ipx][l] = -U[j][i][l];
            if (checkbounds(0, ipx+1, w) && checkbounds(0, ipy, h))
                for (int l = 0; l < 2; l++)
                    V[ipy][ipx+1][l] = -U[j][i][l];
            if (checkbounds(0, ipx, w) && checkbounds(0, ipy+1, h))
                for (int l = 0; l < 2; l++)
                    V[ipy+1][ipx][l] = -U[j][i][l];
            if (checkbounds(0, ipx+1, w) && checkbounds(0, ipy+1, h))
                for (int l = 0; l < 2; l++)
                    V[ipy+1][ipx+1][l] = -U[j][i][l];
        }
}

static void flowinv_iter(float *v, float *u, int w, int h)
{
    float (*V)[w][2] = (void*)v;

    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            float p[2];
            float qx = i + V[j][i][0];
            float qy = j + V[j][i][1];
            bicubic_interpolation(p, u, w, h, 2, qx, qy);
            V[j][i][0] = -p[0];
            V[j][i][1] = -p[1];
        }
}


void adjointflow(
    float *v, // inverse output flow field
    float *u, // input flow field
    int w, // image width
    int h, // image height
    int niter, // number of iterations
    float epsil // tolerance to the boundaries
)
{
    for (int i = 0; i < w*h*2; i++)
        v[i] = -u[i];
    if (epsil > 0)
        flowinv_init(v, u, w, h);

    for (int i = 0; i < niter; i++)
    {
        flowinv_iter(v, u, w, h);
    }
}
