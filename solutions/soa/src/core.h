#ifndef CORE_H
#define CORE_H
// declara las funciones que el main usa
struct soa {
	double x[N];
	double y[N];
	double z[N];
};

void init_pos(struct soa* rxyz, const double rho);
void init_vel(struct soa* vxyz, double* temp, double* ekin);
void forces(const struct soa* rxyz, struct soa* fxyz, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L);
void velocity_verlet(struct soa* rxyz, struct soa* vxyz, struct soa* fxyz, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L);

#endif
