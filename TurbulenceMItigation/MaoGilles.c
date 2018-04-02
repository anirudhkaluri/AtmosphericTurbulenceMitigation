

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "iio.h"
#include "rof.c"
#include "./tvl1flow/mask.c"
#include "./tvl1flow/tvl1flow_lib.c"
#include "./tvl1flow/flowutils.c"
#include "./tvl1flow/viewflow.c"


void temporal_mean(
    float **seq, // Sequence of images
    float *ave, // Temporal average output
    int Nframes, // Number of frames in the sequence
    int w, // frame width
    int h // frame height
)
{

    float *cur=NULL;
    int k,i,j;

    //initialize to 0
    for (i=0; i<w; i++)
    {
        for (j=0; j<h; j++)
        {
            ave[i+j*w]=0.0;
        }
    }

    //sum all the frames
    for (k=0; k<Nframes; k++)
    {
        cur = seq[k];
        for (i=0; i<w; i++)
        {
            for (j=0; j<h; j++)
            {
                ave[i+j*w]=ave[i+j*w]+cur[i+j*w];
            }
        }
    }

    for (i=0; i<w; i++)
    {
        for (j=0; j<h; j++)
        {
            ave[i+j*w]=ave[i+j*w]/Nframes;
        }
    }
    cur = NULL;
    return;
}



static void optflow(
    float *reference, // reference image
    float *frame, // deformed image
    float *vecx, // horizontal optical flow component
    float *vecy, // vertical optical flow component
    int w, // image width
    int h // image height
)
{

    float eps=1e-8; // this variable is used to check how close from 
										// the boundaries of the image we are
    int i,j;
    float pos;

    //set the optical flow parameters (default parameters given by 
		//TVL1opticalflow code)
    float tau = 0.25;
    float lambda = 0.15;
    float theta   = 0.3;
    int nscales = 3;
    float zfactor = 0.5;
    int warps = 1;
    float epsilon = 0.01;
    int verbose = 0;

    //compute the TV-L1 optical flow
    Dual_TVL1_optic_flow_multiscale(frame,reference,vecx,vecy,w,h,tau,lambda,
                                    theta,nscales,zfactor,warps,epsilon,
                                    verbose);

    //we cut vectors which go outside the image domain
    for (i=0; i<w; i++)
    {
        for (j=0; j<h; j++)
        {
            pos=i+vecx[i+j*w];
            if (pos<eps)
            {
                vecx[i+j*w] = eps-i;
            }
            else if (pos>(w-eps))
            {
                vecx[i+j*w] = w-eps-i;
            }

            pos=j+vecy[i+j*w];
            if (pos<eps)
            {
                vecy[i+j*w] = eps-j;
            }
            else if (pos>(h-eps))
            {
                vecy[i+j*w] = h-eps-j;
            }
        }
    }

    return;
}



void MaoGilles_Stabilization(
    float **seq, // input sequence
    float *u, // output restored image
    int width, // image width
    int height, // image height
    int Nframes, // number of frames
    float lambda, // ROF regularization parameter
    float dt, // Bregman time step
		int Nbregman, // number of iterations of the Bregman loop (default 4)
		int Nsplitting // number of iterations of the splitting (default 5)
)
{
    float **seqk=NULL;
    float *seqm=NULL;
    float *du=NULL;
    float *tmp=NULL;
    float *vecx=NULL;
    float *vecy=NULL;
    float **f=NULL;
    float **adjf=NULL;
    int n,i,j,b,iter;
		char name[1024];

    seqk = (float **)xmalloc(sizeof(float *)*Nframes); //Bregman variable
    f = (float **)xmalloc(sizeof(float *)*Nframes); // flow which map u to seq
    adjf = (float **)xmalloc(sizeof(float *)*Nframes); // adjoint flow

    // initialization of the Bregman variable
    for (n=0; n<Nframes; n++)
    {
        seqk[n] = (float *)xmalloc(sizeof(float)*width*height);
        memcpy(seqk[n],seq[n],sizeof(float)*width*height);
    }

    //Memory allocation of flow fields
    vecx = (float *)xmalloc(sizeof(float)*width*height);
    vecy = (float *)xmalloc(sizeof(float)*width*height);
    for (n=0; n<Nframes; n++)
    {
        f[n] = xmalloc(2 * width * height * sizeof(float));
        adjf[n] = xmalloc(2 * width * height * sizeof(float));
    }

    // Auxiliary variable memory allocation
    seqm	= (float *)xmalloc(sizeof(float)*width*height);
    tmp	= (float *)xmalloc(sizeof(float)*width*height);
    du	= (float *)xmalloc(sizeof(float)*width*height);

    // BREGMAN LOOP
    for (b=0; b<Nbregman; b++)
    {
        // Compute the diffeomorphisms between all frames and the reference
        for (n=0; n<Nframes; n++)
        {
            optflow(u,seq[n],vecx,vecy,width,height);
            for (i = 0; i < width * height; i++)
            {
                f[n][2*i] = vecx[i];
                f[n][2*i+1] = vecy[i];
            }
            // compute the adjoint flow
            adjointflow(adjf[n],f[n],width,height,10,0.01);
        }

        //Operator splitting
        for (iter=0; iter<Nsplitting; iter++)
        {
            for (i=0; i<width; i++)
            {
                for (j=0; j<height; j++)
                {
                    du[i+j*width] = 0.0;
                }
            }


            for (n=0; n<Nframes; n++)
            {
                warpwithflow(tmp,f[n],u,width,height,1);// map u on each input
                for (i=0; i<width; i++)
                {
                    for (j=0; j<height; j++)
                    {
                        tmp[i+j*width]=tmp[i+j*width]-seqk[n][i+j*width];
                    }
                }
                warpwithflow(seqm,adjf[n],tmp,width,height,1); //Compute adjoint

                for (i=0; i<width; i++)
                {
                    for (j=0; j<height; j++)
                    {
                        du[i+j*width] = du[i+j*width]+seqm[i+j*width]/Nframes;
                    }
                }
            }

            for (i=0; i<width; i++)
            {
                for (j=0; j<height; j++)
                {
                    du[i+j*width] = u[i+j*width]-dt*du[i+j*width];
                }
            }

            // TV regularization
            RudinOsherFatemi(du,u,width,height,lambda,20);
        } //end splitting loop


        for (n=0; n<Nframes; n++)
        {
            warpwithflow(tmp,f[n],u,width,height,1);

            for (i=0; i<width; i++)
            {
                for (j=0; j<height; j++)
                {
                    seqk[n][i+j*width]=seqk[n][i+j*width]
                        +seq[n][i+j*width]-tmp[i+j*width];
                }
            }
        }
    } //end bregman loop
 
	  // Free memory
    free(du);
    free(tmp);
    free(seqm);
    for (n=0; n<Nframes; n++)
    {
        free(f[n]);
        free(adjf[n]);
        free(seqk[n]);
    }
    free(f);
    free(adjf);
    free(seqk);
    free(vecx);
    free(vecy);

    return;
}
