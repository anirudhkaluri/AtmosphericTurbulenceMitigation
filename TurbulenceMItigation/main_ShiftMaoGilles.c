

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "iio.h"
#include "MaoGilles.c"


int
main(int argc, char **argv)
{

    float **seq=NULL;
    float **seqind=NULL;
    char name[1024];
    float *u=NULL;
		float *vecx=NULL;
    float *vecy=NULL;
    float *f=NULL;
    int Ni,Nf,Nframes,Ns,n,m,i,ind,width,height,Nbregman,Nsplitting,Saveflow,save;
    float lambda;
    float dt;

    time_t start,end;
    float t;


    if (argc<10) {
        fprintf(stderr,"usage: ShiftMaoGilles image_generic_name first_index \
                last_index lambda delta temporal_support_length \
								saveflow_0_or_1\n\n");
        exit(EXIT_FAILURE);
    }

    // Get the input parameters
    Ni = atoi(argv[2]); //get the index of the first frame to load
    Nf = atoi(argv[3]); //get the index of the last frame to load
    Nframes = Nf-Ni+1; //compute the number of frames in the sequence
    lambda = atof(argv[4]); //regularization parameter
    dt = atof(argv[5]); //time step
    Ns = atoi(argv[6]); //temporal support length
		Nbregman = atoi(argv[7]); //number of iteration for the Bregman loop
		Nsplitting = atoi(argv[8]); //number of iteration for the Splitting loop
		Saveflow = atoi(argv[9]); //1 of you want to save the flows between the 
															//first restored frame and the first ten input 
															//frames, 0 otherwise

    fprintf(stderr,"Ni=%d - Nf=%d - Ns = %d - Lambda=%f - dt=%f - Nbregman=%d \
						- Nsplitting=%d - Saveflow=%d\n",Ni,Nf,Ns,lambda,dt,Nbregman, \
						Nsplitting,Saveflow);


    //Input Sequence Memory Allocations
    seq	= (float **)xmalloc(sizeof(float *)*Nframes);

    //load sequence
    for (n=Ni; n<Nf+1; n++)
    {
        sprintf(name,argv[1],n);
        seq[n-Ni] = iio_read_image_float(name,&width,&height);
    }
    fprintf(stderr,"Sequence loaded\n");

		uint8_t (**view)[3] = matrix_build(width, height, sizeof**view);

#pragma omp parallel private(seqind,u,vecx,vecy,f,n,m,i,ind,name)
    {

#pragma omp for
        for (n=0;n<(Nframes-Ns+1);n++)
        {
				    u = (float *)xmalloc(sizeof(float)*width*height); //Output image
    				seqind	= (float **)xmalloc(sizeof(float *)*Ns); //temporal window sequence

            fprintf(stderr,"\n *** Window %d ***\n",n);
            // Assign the correct input images to the current temporal window
            for (ind=0;ind<Ns;ind++) seqind[ind] = seq[n+ind];

						//initialization of u by the temporal mean
						temporal_mean(seqind,u,Ns,width,height);

            // Perform the restoration process
						MaoGilles_Stabilization(seqind,u,width,height,Ns,lambda,dt,\
																		Nbregman,Nsplitting);

            // Save the restored image
            sprintf(name,"Restored_%03d.png",n+117);
            iio_save_image_float(name, u, width, height);

						// Save the optical flow between the first restored image and the
						// ten first input frames
						if ((Saveflow==1) && (n==0)) {
							vecx = (float *)xmalloc(sizeof(float)*width*height);
							vecy = (float *)xmalloc(sizeof(float)*width*height);
					 		f = xmalloc(2 * width * height * sizeof(float));

							for (m=0; m<(Nframes>10?10:Nframes); m++)
        			{
            		optflow(u,seq[m],vecx,vecy,width,height);
            		for (i = 0; i < width * height; i++)
            		{
                	f[2*i] = vecx[i];
                	f[2*i+1] = vecy[i];
            		}

								sprintf(name,"flow%03d.png",m+1);
								viewflow_flat(view[0][0], f, width, height, 1);
								iio_save_image_uint8_vec(name, view[0][0], width, height, 3);
        			}
   						free(vecx);
							free(vecy);
							free(f);
						}
				   free(u);
				   free(seqind);
        }
    }

    fprintf(stderr,"Sequence saved\n");

    // Memory freeing
    free(seq);
    exit(EXIT_SUCCESS);
}
