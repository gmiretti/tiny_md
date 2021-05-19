#include "parameters.h"

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6))) 

double minimum_image_ispc(double cordi, const double cell_length)
{
    // imagen más cercana-esto le dice estilo que si la caja de sim. va de 0 a cell_length, 
    // con este cálculo veo qué partículas considera como dentro de la celda de una partíc. (ver 1.1.2 TinyMD.pdf)

    if (cordi <= -0.5 * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5 * cell_length) {
        // printf("minimum_image(%f, %f)\n", cordi, cell_length); // for debug
        cordi -= cell_length;
    }
    return cordi;
    // return (cordi - ((cordi > 0)?1:-1) * ((abs(cordi)>0.5*cell_length)?cell_length:0)); // worse
}



export void forces_ispc(uniform double rxyz[], uniform double fxyz[], 
            uniform double* uniform epot, uniform double* uniform pres, uniform double* uniform temp, 
            uniform const double rho, uniform const double V, uniform const double L)
{
    // calcula las fuerzas LJ (12-6)

    foreach (i = 0 ... 3 * N) {
        fxyz[i] = 0.0;
    }
    *epot = 0.0;
    uniform double pres_vir = 0.0;
    uniform double rcut2 = RCUT * RCUT;

    foreach (jdx = 0 ... N, idx = 0 ... N - 1) {
        varying int i = idx *3;
        varying int j = jdx * 3;

        varying double xi = rxyz[i + 0];
        varying double yi = rxyz[i + 1];
        varying double zi = rxyz[i + 2];

        if (j >= i + 3) {

            varying double xj = rxyz[j + 0];
            varying double yj = rxyz[j + 1];
            varying double zj = rxyz[j + 2];

            // distancia mínima entre r_i y r_j
            varying double rx = xi - xj;
            // printf("minimum_image(%f, %f) =", rx, L); // for debug
            rx = minimum_image_ispc(rx, L);
            // printf("%f \n", rx); // for debug
            varying double ry = yi - yj;
            ry = minimum_image_ispc(ry, L);
            varying double rz = zi - zj;
            rz = minimum_image_ispc(rz, L);

            varying double rij2 = rx * rx + ry * ry + rz * rz;

            // con esto se fija si está dentro del radio de corte. (con módulo de cada vector posic.)
            if (rij2 <= rcut2) {
                varying double r2inv = 1.0 / rij2;
                varying double r6inv = r2inv * r2inv * r2inv;

                varying double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

                fxyz[i + 0] += fr * rx;
                fxyz[i + 1] += fr * ry;
                fxyz[i + 2] += fr * rz;

                fxyz[j + 0] -= fr * rx;
                fxyz[j + 1] -= fr * ry;
                fxyz[j + 2] -= fr * rz;

                *epot = reduce_add(4.0 * r6inv * (r6inv - 1.0) - ECUT);
                pres_vir = reduce_add(fr * rij2);
            }
        }
    }
    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;
}
