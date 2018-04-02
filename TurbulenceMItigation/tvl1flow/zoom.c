


#ifndef ZOOM_H
#define ZOOM_H

#include <assert.h>


void zoom_size(
    int nx,             //width of the orignal image
    int ny,             //height of the orignal image
    int *nxx,           //width of the zoomed image
    int *nyy,           //height of the zoomed image
    float factor  //zoom factor between 0 and 1
)
{
    *nxx = (int)((float) nx * factor + 0.5);
    *nyy = (int)((float) ny * factor + 0.5);
}



void zoom(
    const float *I,    //input image
    float *Iout,       //output image
    int nx,            //image width
    int ny,            //image height
    float factor //zoom factor between 0 and 1
)
{
    int nxx, nyy;

    zoom_size(nx, ny, &nxx, &nyy, factor);

    for (int i = 0; i < nyy; i++)

        for (int j = 0; j < nxx; j++)
        {
            const int l =(int)((float)i/factor);
            const int k =(int)((float)j/factor);

            int mx = k - 1;
            int dx = k + 1;
            int my = l - 1;
            int dy = l + 1;

            if(mx < 0) mx = 0;
            if(dx >= nx) dx = nx - 1;
            if(my < 0) my = 0;
            if(dy >= ny) dy = ny - 1;

            //smooth with the neighbors
            float value =
                0.07842776544 * (I[mx + nx * my] + I[mx + nx * dy]
                                    + I[dx + nx * my] + I[dx + nx * dy]) +
                0.1231940459  * (I[k + nx * my]  + I[mx + nx * l]
                                    + I[dx + nx * l]  + I[k + nx * dy])  +
                0.1935127547  *  I[k + nx * l];

                Iout[i * nxx + j] = value;
        }
}



void zoom_out(
    const float *I, //input image
    float *Iout,    //output image
    int nx,         //width of the original image
    int ny,         //height of the original image
    int nxx,        //width of the zoomed image
    int nyy         //height of the zoomed image
)
{
    const float factorx = ((float)nxx / nx);
    const float factory = ((float)nyy / ny);

    for (int i = 0; i < nyy; i++)

        for (int j = 0; j < nxx; j++)
        {
            int l =  (int) ((float)i/factory);
            int k =  (int) ((float)j/factorx);

            assert(l < ny);
            assert(k < nx);
            Iout[i * nxx + j] = I[k + nx * l];
        }
}



#endif
