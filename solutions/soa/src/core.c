#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6))) 


void init_pos(struct soa* rxyz, const double rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0)); //cbrt es cubic root, ceil(x) vuelve int a x
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
		// 4 calculos de ryxz porque son 4 direcc. de vectores red para FCC
		rxyz-> x[idx] = i * a; // x
		rxyz-> y[idx] = j * a; // y 
		rxyz-> z[idx] = k * a; // z
                // rxyz[idx + 0] = i * a; // x
                // rxyz[idx + 1] = j * a; // y
                // rxyz[idx + 2] = k * a; // z
                    // del mismo átomo
		rxyz-> x[idx + 1] = (i + 0.5) * a;
		rxyz-> y[idx + 1] = (j + 0.5) * a;
		rxyz-> z[idx + 1] = k * a;
                // rxyz[idx + 3] = (i + 0.5) * a;
                // rxyz[idx + 4] = (j + 0.5) * a;
                // rxyz[idx + 5] = k * a;

		rxyz-> x[idx + 2] = (i + 0.5) * a;
		rxyz-> y[idx + 2] = j * a;
		rxyz-> z[idx + 2] = (k + 0.5) * a;
                // rxyz[idx + 6] = (i + 0.5) * a;
                // rxyz[idx + 7] = j * a;
                // rxyz[idx + 8] = (k + 0.5) * a;

		rxyz-> x[idx + 3] = i * a;
		rxyz-> y[idx + 3] = (j + 0.5) * a;
		rxyz-> z[idx + 3] = (k + 0.5) * a;
                // rxyz[idx + 9] = i * a;
                // rxyz[idx + 10] = (j + 0.5) * a;
                // rxyz[idx + 11] = (k + 0.5) * a;

                idx += 4; // 12
            }
        }
    }
}


void init_vel(struct soa* vxyz, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    // Struct. of arrays
    for (int i = 0; i < N; i++) {
        vxyz-> x[i] = rand() / (double)RAND_MAX - 0.5;
        vxyz-> y[i] = rand() / (double)RAND_MAX - 0.5;
        vxyz-> z[i] = rand() / (double)RAND_MAX - 0.5;

        sumvx += vxyz-> x[i];
        sumvy += vxyz-> y[i];
        sumvz += vxyz-> z[i];
        sumv2 += vxyz-> x[i] * vxyz-> x[i] + vxyz-> y[i] * vxyz-> y[i]
            + vxyz-> z[i] * vxyz-> z[i];
    }


    // un solo arreglo con todas las velocidades
    // for (int i = 0; i < 3 * N; i += 3) {
    //     vxyz[i + 0] = rand() / (double)RAND_MAX - 0.5;
    //     vxyz[i + 1] = rand() / (double)RAND_MAX - 0.5;
    //     vxyz[i + 2] = rand() / (double)RAND_MAX - 0.5;

    //     sumvx += vxyz[i + 0];
    //     sumvy += vxyz[i + 1];
    //     sumvz += vxyz[i + 2];
    //     sumv2 += vxyz[i + 0] * vxyz[i + 0] + vxyz[i + 1] * vxyz[i + 1]
    //         + vxyz[i + 2] * vxyz[i + 2];
    // }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);
    
    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        vxyz-> x[i] = (vxyz-> x[i] - sumvx) * sf;
        vxyz-> y[i] = (vxyz-> y[i] - sumvy) * sf;
        vxyz-> z[i] = (vxyz-> z[i] - sumvz) * sf;
    }

    // for (int i = 0; i < 3 * N; i += 3) { // elimina la velocidad del centro de masa
    //     // y ajusta la temperatura
    //     vxyz[i + 0] = (vxyz[i + 0] - sumvx) * sf;
    //     vxyz[i + 1] = (vxyz[i + 1] - sumvy) * sf;
    //     vxyz[i + 2] = (vxyz[i + 2] - sumvz) * sf;
    // }
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


void forces(const struct soa* rxyz, struct soa* fxyz, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        fxyz-> x[i] = 0.0;
	fxyz-> y[i] = 0.0;
	fxyz-> z[i] = 0.0;
    }
    // for (int i = 0; i < 3 * N; i++) {
    //     fxyz[i] = 0.0;
    // }
    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < (N - 1); i++) {
    // for (int i = 0; i < 3 * (N - 1); i += 3) {

        double xi = rxyz-> x[i];
        double yi = rxyz-> y[i];
        double zi = rxyz-> z[i];

        for (int j = i + 1; j < N; j++) {

            double xj = rxyz-> x[j];
            double yj = rxyz-> y[j];
            double zj = rxyz-> z[j];

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

                fxyz-> x[i] += fr * rx;
                fxyz-> y[i] += fr * ry;
                fxyz-> z[i] += fr * rz;

                fxyz-> x[j] -= fr * rx;
                fxyz-> y[j] -= fr * ry;
                fxyz-> z[j] -= fr * rz;

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


void velocity_verlet(struct soa* rxyz, struct soa* vxyz, struct soa* fxyz, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        rxyz-> x[i] += vxyz-> x[i] * DT + 0.5 * fxyz-> x[i] * DT * DT;
        rxyz-> y[i] += vxyz-> y[i] * DT + 0.5 * fxyz-> y[i] * DT * DT;
        rxyz-> z[i] += vxyz-> z[i] * DT + 0.5 * fxyz-> z[i] * DT * DT;

        rxyz-> x[i] = pbc(rxyz-> x[i], L);
        rxyz-> y[i] = pbc(rxyz-> y[i], L);
        rxyz-> z[i] = pbc(rxyz-> z[i], L);

        vxyz-> x[i] += 0.5 * fxyz-> x[i] * DT;
        vxyz-> y[i] += 0.5 * fxyz-> y[i] * DT;
        vxyz-> z[i] += 0.5 * fxyz-> z[i] * DT;
    }

    forces(rxyz, fxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        vxyz-> x[i] += 0.5 * fxyz-> x[i] * DT;
        vxyz-> y[i] += 0.5 * fxyz-> y[i] * DT;
        vxyz-> z[i] += 0.5 * fxyz-> z[i] * DT;

        sumv2 += vxyz-> x[i] * vxyz-> x[i] + vxyz-> y[i] * vxyz-> y[i]
            + vxyz-> z[i] * vxyz-> z[i];
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
