#define _XOPEN_SOURCE 500  // M_PI
#include "core.h"
#include "parameters.h"
#include "wtime.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main()
{

    // archivos y def. variables
    FILE *file_xyz, *file_thermo;
    file_xyz = fopen("trajectory.xyz", "w");
    file_thermo = fopen("thermo.log", "w");
    double Ekin, Epot, Temp, Pres; // variables macroscopicas
    double Rho, cell_V, cell_L, tail, Etail, Ptail;
    // double *rxyz, *vxyz, *fxyz; // variables microscopicas
   
    // Structure of arrays-útil para vectorizar y operar de a muchos elementos. 
    // Además esto suele ayudar al autovectorizador.
    // ocupa el tam. de un puntero (8 bytes-esto si usa lo de soa *rxyz)
    // + el tam. de calloc en el heat/p que es el tam. de los struct... N doubles por coord.+8 bytes entonces.
    // tonces: todas las opciones ocupan lo mismo porque total es un solo puntero
    // struct soa {
    //         double x[N];
    //         double y[N];
    //         double z[N];
    // };

    struct soa *rxyz = calloc(1, sizeof(struct soa));
    struct soa *vxyz = calloc(1, sizeof(struct soa));
    struct soa *fxyz = calloc(1, sizeof(struct soa));
    // struct soa rxyz; // lo pone en el ¿heap/t? con struct soa *rxyz = calloc(1, sizeof(struct soa));
    // rxyz.x[i]; // rxyz-> x[i]

    // pide memoria-no se si significa esto tanto para vers. original comentada y soa tb. :P
    // rxyz = (double*)malloc(3 * N * sizeof(double));
    // vxyz = (double*)malloc(3 * N * sizeof(double));
    // fxyz = (double*)malloc(3 * N * sizeof(double));

    // imprime en pantalla-tal vez comentar
    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);
    printf("# densidad, volumen, energía potencial media, presión media\n");
    fprintf(file_thermo, "# t Temp Pres Epot Etot\n"); // el encabezado del archivo de salida thermo.log

    srand(SEED); // dynamic initialization (de SEED)-266-ModernC.pdf
    double t = 0.0, sf;
    double Rhob;
    Rho = RHOI;				// RHOI es densidad inicial
    init_pos(rxyz, Rho);		// inicializa posiciones
    double start = wtime();		// para llevar cuenta del tiempo de corrida
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1 * (double)m;
        cell_V = (double)N / Rho;	// volumen celda
        cell_L = cbrt(cell_V);		// largo celda
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0; // para corregir cálculo de E y P
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);
        for (int k = 0; k < N; k++) { // reescaleo posiciones a nueva densidad
            rxyz-> x[k] *= sf;
            rxyz-> y[k] *= sf;
            rxyz-> z[k] *= sf;
        }

        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < TEQ; i++) { // loop de equilibracion-aca no guarda propiedades termodin.
	    // actualiza posiciones y propiedades con Velocity Verlet
            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                vxyz-> x[k] *= sf;
		vxyz-> y[k] *= sf;
		vxyz-> z[k] *= sf;
            }
        }

        int mes = 0;
        double epotm = 0.0, presm = 0.0;
        for (i = TEQ; i < TRUN; i++) { // loop de medicion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
               vxyz-> x[k] *= sf;
	       vxyz-> y[k] *= sf;
	       vxyz-> z[k] *= sf;
            }           

            if (i % TMES == 0) {	// escritura de archivos
                Epot += Etail;		// suma factores para ir corrigiendo la E y otras calculadas
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;
		
		
                fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot, Epot + Ekin);
                fprintf(file_xyz, "%d\n\n", N);
                for (int k = 0; k < N; k++) {
                    fprintf(file_xyz, "Ar %e %e %e\n", rxyz-> x[k], rxyz-> y[k], rxyz-> z[k]);
                }
            }

            t += DT;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (double)mes, presm / (double)mes);
    }
    // poner T como parámetro

    double elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    return 0;
}
