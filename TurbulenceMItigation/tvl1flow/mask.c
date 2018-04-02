

#ifndef MASK_H
#define MASK_H


void divergence(
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

        //apply the divergence to the first and last rows
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
 * (see [2] for details)
 *
 **/
void forward_gradient(
    const float *f, //input image
    float *fx,      //computed x derivative
    float *fy,      //computed y derivative
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
            }
        }

        //apply the gradient to the last row
        for(int j = 0; j < nx-1; j++)
        {
            const int p = (ny-1) * nx + j;

            fx[p] = f[p+1] - f[p];
            fy[p] = 0;
        }

        //apply the gradient to the last column
        for(int i = 1; i < ny; i++)
        {
            const int p = i * nx-1;

            fx[p] = 0;
            fy[p] = f[p+nx] - f[p];
        }

        fx[ny * nx - 1] = 0;
        fy[ny * nx - 1] = 0;
}


/**
 *
 * Function to compute the gradient with centered differences
 *
 **/
void centered_gradient(
    const float *input,  //input image
    float *dx,           //computed x derivative
    float *dy,           //computed y derivative
    const int nx,        //image width
    const int ny         //image height
)
{
        //apply the gradient to the center body of the image
        for(int i = 1; i < ny-1; i++)
        {
            for(int j = 1; j < nx-1; j++)
            {
                const int k = i * nx + j;
                dx[k] = 0.5*(input[k+1] - input[k-1]);
                dy[k] = 0.5*(input[k+nx] - input[k-nx]);
            }
        }

        //apply the gradient to the first and last rows
        for(int j = 1; j < nx-1; j++)
        {
            dx[j] = 0.5*(input[j+1] - input[j-1]);
            dy[j] = 0.5*(input[j+nx] - input[j]);

            const int k = (ny - 1) * nx + j;

            dx[k] = 0.5*(input[k+1] - input[k-1]);
            dy[k] = 0.5*(input[k] - input[k-nx]);
        }

        //apply the gradient to the first and last columns
        for(int i = 1; i < ny-1; i++)
        {
            const int p = i * nx;
            dx[p] = 0.5*(input[p+1] - input[p]);
            dy[p] = 0.5*(input[p+nx] - input[p-nx]);

            const int k = (i+1) * nx - 1;

            dx[k] = 0.5*(input[k] - input[k-1]);
            dy[k] = 0.5*(input[k+nx] - input[k-nx]);
        }

        //apply the gradient to the four corners
        dx[0] = 0.5*(input[1] - input[0]);
        dy[0] = 0.5*(input[nx] - input[0]);

        dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
        dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

        dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
        dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

        dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
        dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}
#endif
