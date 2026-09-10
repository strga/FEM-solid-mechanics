#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "csr.h"
#include "mesh.h"

// Configuration structure for simulation parameters
typedef struct {
    char simulationParametersFile[256];
    char displacementFile[256];
    char energyFile[256];
    char energyDampingFile[256];
    char timeSeriesFile[256];
    char vtkDirectory[256];
} SimulationConfig;

// Function prototypes
void initializeSimulationConfig(SimulationConfig *config);
int findNeumannBoundaryNode(const mesh *M);

double calculateWorkByReference(double *u, double *U, csr *K, int nDOF);
double calculateInitialEnergy(double *u_s, csr *K, int nDOF);
// double calculateInitialEnergy(const double *u0, const double *v0, const csr *K, const csr *M, int nDOF);
double calculateKineticEnergy(double *v, csr *M, int nDOF);
double calculatePotentialEnergy(double *u, csr *K, int nDOF);

void NewmarkIntegrate(mesh *this, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, double dt, int nSteps, int useDamping, int saveVibration, int saveData, SimulationConfig *config);
void NewmarkIntegrate_Nonlinear(mesh *M, csr *M_csr, csr *D_csr, double *F_ext, double *u, double *v, double *a, double dt, int nSteps, int maxNewton, double tol, int useDamping, int saveVibration, int saveData, SimulationConfig *config, double (*psi_force_3D)(double,double,double,double,double,double,double,double,double,double*,double*,double*), double *dirichletValues, double F_mag, double azimuth_deg, double elevation_deg);

void Mesh_Discretize_NonlinearDynamics3D(mesh *M, triplet *stiff, double *residual, double *f_ext, double *P_int, double *u_global, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg);
#endif // DYNAMICS_H
