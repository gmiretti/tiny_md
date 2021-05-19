#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6))) 


void init_pos(struct aos* rxyz, const double rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0)); //cbrt es cubic root, ceil(x) vuelve int a x
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
		// 4 calculos de ryxz porque son 4 direcc. de vectores red para FCC
                rxyz[idx].x = i * a; // x
                rxyz[idx].y = j * a; // y
                rxyz[idx].z = k * a; // z
                // rxyz[idx + 0] = i * a; // x
                // rxyz[idx + 1] = j * a; // y
                // rxyz[idx + 2] = k * a; // z
                    // del mismo átomo
                rxyz[idx + 1].x = (i + 0.5) * a;
                rxyz[idx + 1].y = (j + 0.5) * a;
                rxyz[idx + 1].z = k * a;
                // rxyz[idx + 3] = (i + 0.5) * a;
                // rxyz[idx + 4] = (j + 0.5) * a;
                // rxyz[idx + 5] = k * a;

                rxyz[idx + 2].x = (i + 0.5) * a;
                rxyz[idx + 2].y = j * a;
                rxyz[idx + 2].z = (k + 0.5) * a;
                // rxyz[idx + 6] = (i + 0.5) * a;
                // rxyz[idx + 7] = j * a;
                // rxyz[idx + 8] = (k + 0.5) * a;

                rxyz[idx + 3].x = i * a;
                rxyz[idx + 3].y = (j + 0.5) * a;
                rxyz[idx + 3].z = (k + 0.5) * a;
                // rxyz[idx + 9] = i * a;
                // rxyz[idx + 10] = (j + 0.5) * a;
                // rxyz[idx + 11] = (k + 0.5) * a;

                idx += 4;
            }
        }
    }
}


void init_vel(struct aos* vxyz, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    // un solo arreglo con todas las velocidades
    for (int i = 0; i < N; i++) {
        vxyz[i].x = rand() / (double)RAND_MAX - 0.5;
        vxyz[i].y = rand() / (double)RAND_MAX - 0.5;
        vxyz[i].z = rand() / (double)RAND_MAX - 0.5;

        sumvx += vxyz[i].x;
        sumvy += vxyz[i].y;
        sumvz += vxyz[i].z;
        sumv2 += vxyz[i].x * vxyz[i].x + vxyz[i].y * vxyz[i].y
            + vxyz[i].z * vxyz[i].z;
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        vxyz[i].x = (vxyz[i].x - sumvx) * sf;
        vxyz[i].y = (vxyz[i].y - sumvy) * sf;
        vxyz[i].z = (vxyz[i].z - sumvz) * sf;
    }
}


static double minimum_image(double cordi, const double cell_length)
{
    // imagen más cercana-esto le dice estilo que si la caja de sim. va de 0 a cell_length, 
    // con este cálculo veo qué partículas considera como dentro de la celda de una partíc. (ver 1.1.2 TinyMD.pdf)

    if (cordi <= -0.5 * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5 * cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}


void forces(const struct aos* rxyz, struct aos* fxyz, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        fxyz[i].x = 0.0;
        fxyz[i].y = 0.0;
        fxyz[i].z = 0.0;
    }
    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < (N - 1); i++) {

        double xi = rxyz[i].x;
        double yi = rxyz[i].y;
        double zi = rxyz[i].z;

        for (int j = i + 1; j < N; j++) {

            double xj = rxyz[j].x;
            double yj = rxyz[j].y;
            double zj = rxyz[j].z;

            // distancia mínima entre r_i y r_j
            double rx = xi - xj;
            rx = minimum_image(rx, L);
            double ry = yi - yj;
            ry = minimum_image(ry, L);
            double rz = zi - zj;
            rz = minimum_image(rz, L);

            double rij2 = rx * rx + ry * ry + rz * rz;

	    // con esto se fija si está dentro del radio de corte. (con módulo de cada vector posic.)
            if (rij2 <= rcut2) {
                double r2inv = 1.0 / rij2;
                double r6inv = r2inv * r2inv * r2inv;

                double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

                fxyz[i].x += fr * rx;
                fxyz[i].y += fr * ry;
                fxyz[i].z += fr * rz;

                fxyz[j].x -= fr * rx;
                fxyz[j].y -= fr * ry;
                fxyz[j].z -= fr * rz;

                *epot += 4.0 * r6inv * (r6inv - 1.0) - ECUT;
                pres_vir += fr * rij2;
            }
        }
    }
    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;
}


static double pbc(double cordi, const double cell_length)
{
    // condiciones periodicas de contorno coordenadas entre [0,L)
    if (cordi <= 0) {
        cordi += cell_length;
    } else if (cordi > cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}


void velocity_verlet(struct aos* rxyz, struct aos* vxyz, struct aos* fxyz, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        rxyz[i].x += vxyz[i].x * DT + 0.5 * fxyz[i].x * DT * DT;
        rxyz[i].y += vxyz[i].y * DT + 0.5 * fxyz[i].y * DT * DT;
        rxyz[i].z += vxyz[i].z * DT + 0.5 * fxyz[i].z * DT * DT;

        rxyz[i].x = pbc(rxyz[i].x, L);
        rxyz[i].y = pbc(rxyz[i].y, L);
        rxyz[i].z = pbc(rxyz[i].z, L);

        vxyz[i].x += 0.5 * fxyz[i].x * DT;
        vxyz[i].y += 0.5 * fxyz[i].y * DT;
        vxyz[i].z += 0.5 * fxyz[i].z * DT;
    }

    forces(rxyz, fxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        vxyz[i].x += 0.5 * fxyz[i].x * DT;
        vxyz[i].y += 0.5 * fxyz[i].y * DT;
        vxyz[i].z += 0.5 * fxyz[i].z * DT;

        sumv2 += vxyz[i].x * vxyz[i].x + vxyz[i].y * vxyz[i].y
            + vxyz[i].z * vxyz[i].z;
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
