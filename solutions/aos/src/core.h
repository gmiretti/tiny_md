#ifndef CORE_H
#define CORE_H
// declara las funciones que el main usa
struct aos {
	double x, y, z;
};

void init_pos(struct aos* rxyz, const double rho);
void init_vel(struct aos* vxyz, double* temp, double* ekin);
void forces(const struct aos* rxyz, struct aos* fxyz, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L);
void velocity_verlet(struct aos* rxyz, struct aos* vxyz, struct aos* fxyz, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L);

#endif
