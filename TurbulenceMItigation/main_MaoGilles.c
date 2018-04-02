
#include <stdio.h>
#include <stdlib.h>
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
    char name[1024];
    float *u=NULL;
		float *umean=NULL;
		float *vecx=NULL;
    float *vecy=NULL;
    float *f=NULL;
    int Ni,Nf,Nframes,n,m,i,width,height,Nbregman,Nsplitting,Saveflow;
    float lambda;
    float dt;
		
		clock_t start;
		clock_t end;


    if (argc<9) {
        fprintf(stderr,"usage: MaoGilles image_generic_name first_index \
                last_index lambda delta Nbregman Nsplitting \
								saveflow_0_or_1\n\n");
        exit(EXIT_FAILURE);
    }

    // Get the input parameters
    Ni = atoi(argv[2]); //get the index of the first frame to load
    Nf = atoi(argv[3]); //get the index of the last frame to load
    Nframes = Nf-Ni+1; //compute the number of frames in the sequence
    lambda = atof(argv[4]); //regularization parameter
    dt = atof(argv[5]); //time step
		Nbregman = atoi(argv[6]); //number of iteration for the Bregman loop
		Nsplitting = atoi(argv[7]); //number of iteration for the Splitting loop
		Saveflow = atoi(argv[8]); //1 of you want to save the flows between the 
															//first restored frame and the first ten input 
															//frames, 0 otherwise

    fprintf(stderr,"Ni=%d - Nf=%d - Lambda=%f - dt=%f - Nbregman=%d - \
						Nsplitting=%d - Saveflow=%d\n",Ni,Nf,lambda,dt,Nbregman,\
						Nsplitting,Saveflow);

    //Memory Allocations
    seq	= (float **)xmalloc(sizeof(float *)*Nframes); //Original sequence

    //load sequence
    for (n=Ni; n<Nf+1; n++)
    {
        sprintf(name,argv[1],n);
        seq[n-Ni] = iio_read_image_float(name,&width,&height);
    }
    fprintf(stderr,"Sequence loaded\n");

    // Allocation of the output image
    u = (float *)xmalloc(sizeof(float)*width*height);
		umean = (float *)xmalloc(sizeof(float)*width*height);

    //Initialize u with the temporal mean
    temporal_mean(seq,u,Nframes,width,height);
		iio_save_image_float("Mean.png", u, width, height);

    // Perform the restoration process
		start=clock();
    MaoGilles_Stabilization(seq,u,width,height,Nframes,lambda,dt,Nbregman,\
														Nsplitting);
		end=clock();

    // Save the restored image
    iio_save_image_float("Restored.png", u, width, height);
    fprintf(stderr,"Image saved\n");
		fprintf(stderr,"Computed in %ld s\n",(end-start)/ CLOCKS_PER_SEC);

		// Save the optical flow between the first restored image and the
		// ten first input frames
		if (Saveflow==1) {
			uint8_t (**view)[3] = matrix_build(width, height, sizeof**view);
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
    // Memory freeing
    free(u);
    free(seq);

    exit(EXIT_SUCCESS);
}
