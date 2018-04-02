

#ifndef WARP_H
#define WARP_H


#define BOUNDARY_CONDITION 0


int neumann_bc(int x, int nx, float *mask)
{

    if(x < 0)
    {
        x = 0;
        *mask = 0;
    }
    else if (x >= nx)
    {
        x = nx - 1;
        *mask = 0;
    }

    return x;
}

int periodic_bc(int x, int nx, float *mask)
{
    if(x < 0)
    {
        const int n   = 1 - (int)(x/(nx+1));
        const int ixx = x + n * nx;

        x =   ixx% nx;
        *mask = 0;
    }
    else if(x >= nx)
    {
        x = x % nx;
        *mask = 0;
    }

    return x;
}


int symmetric_bc(int x, int nx, float *mask)
{
    if(x < 0)
    {
        const int borde = nx - 1;
        const int xx = -x;
        const int n  = (int)(xx/borde) % 2;

        if ( n ) x = borde - ( xx % borde );
        else x = xx % borde;
        *mask = 0;
    }
    else if ( x >= nx )
    {
        const int borde = nx - 1;
        const int n = (int)(x/borde) % 2;

        if ( n ) x = borde - ( x % borde );
        else x = x % borde;
        *mask = 0;
    }

    return x;
}

/**
 *
 * Function to warp the image using bilinear interpolation
 *
 **/
void bilinear_interpolation(
    const float *input, //image to be warped
    const float *u,     //x component of the vector field
    const float *v,     //y component of the vector field
    float *output,      //warped output image
    const int nx,       //width of the image
    const int ny,       //height of the image
    float *mask         //mask to detect the motions outside the image
)
{
    for(int i = 0; i < nx * ny; i++)
        mask[i] = 1;

        for(int j = 0; j < ny; j++)

            for(int i = 0; i < nx; i++)
            {

                const float uu = (float) (i + u[i + nx * j]);
                const float vv = (float) (j + v[i + nx * j]);

                const int sx = (uu < 0)? -1: 1;
                const int sy = (vv < 0)? -1: 1;
                const int p  = i + nx * j;

                int x, y, dx, dy;

                switch(BOUNDARY_CONDITION)
                {

                case 0:
                    x  = neumann_bc((int) uu, nx, mask+p);
                    y  = neumann_bc((int) vv, ny, mask+p);
                    dx = neumann_bc(x + sx, nx,   mask+p);
                    dy = neumann_bc(y + sy, ny,   mask+p);
                    break;

                case 1:
                    x  = periodic_bc((int) uu, nx, mask+p);
                    y  = periodic_bc((int) vv, ny, mask+p);
                    dx = periodic_bc(x + sx, nx,   mask+p);
                    dy = periodic_bc(y + sy, ny,   mask+p);
                    break;

                case 2:
                    x  = symmetric_bc((int) uu, nx, mask+p);
                    y  = symmetric_bc((int) vv, ny, mask+p);
                    dx = symmetric_bc(x + sx, nx,   mask+p);
                    dy = symmetric_bc(y + sy, ny,   mask+p);
                    break;

                default:
                    x  = neumann_bc((int) uu, nx, mask+p);
                    y  = neumann_bc((int) vv, ny, mask+p);
                    dx = neumann_bc(x + sx, nx,   mask+p);
                    dy = neumann_bc(y + sy, ny,   mask+p);
                    break;
                }

                const float p1 = input[x  + nx * y];
                const float p2 = input[dx + nx * y];
                const float p3 = input[x  + nx * dy];
                const float p4 = input[dx + nx * dy];

                const float e1 = ((float) sx * (uu - x));
                const float E1 = ((float) 1.0 - e1);
                const float e2 = ((float) sy * (vv - y));
                const float E2 = ((float) 1.0 - e2);

                const float w1 = E1 * p1 + e1 * p2;
                const float w2 = E1 * p3 + e1 * p4;

                output[i + nx * j] = E2 * w1 + e2 * w2;
            }
}


#endif
