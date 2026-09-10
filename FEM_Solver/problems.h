#ifndef _PROBLEMS_H
#define _PROBLEMS_H

#include "mesh.h"
#include "triplet.h"
#include "csr.h"
#include "dynamics.h"

// Function declarations
void InitializePoissonProblem(mesh *M, triplet *stiff, csr *A_csr, double **b, double **dirichletValues, double **x, IntegrationType integrationType);
void ExecutePoissonAnalysis(mesh *M, csr *A_csr, double *b, double *x, int nDOF, int solverType, double tolerance);
void CleanupPoissonProblem(triplet *stiff, csr *A_csr, double *b, double *dirichletValues, double *x);
/* ------------------------------------------------------------------------------------------------- */
void InitializePoissonProblem3D(mesh *M, triplet *stiff, csr *A_csr, double **b, double **dirichletValues, double **x);
void ExecutePoissonAnalysis3D(mesh *M, csr *A_csr, double *b, double *x, int nDOF, int solverType, double tolerance);
void CleanupPoissonProblem3D(triplet *stiff, csr *A_csr, double *b, double *dirichletValues, double *x);
/* ------------------------------------------------------------------------------------------------- */
void InitializeStaticsProblem(mesh *M, triplet *stiff, csr *K_csr, double **b, double **dirichletValues, double **x, double angle_degrees, int gravity);
void ExecuteStaticsAnalysis(mesh *M, csr *K_csr, double *b, double *x, int nDOF, int solverType, double tolerance);
void CleanupStaticsProblem(triplet *stiff, csr *K_csr, double *b, double *dirichletValues, double *x);
/* ------------------------------------------------------------------------------------------------- */
void InitializeDynamicProblem(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, int *nDOF, double alpha, double eps1, double eps2, int gravity, SimulationConfig *config);
void ExecuteDynamicAnalysis(mesh *M, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, int nDOF, double omega, double eps1, double eps2, double dt, int nSteps, int Damping, int saveVKT, int saveData, SimulationConfig *config);
void CleanupDynamicProblem(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double *u, double *v, double *a, double *F, double *dirichletValues);
/* ------------------------------------------------------------------------------------------------- */
void InitializeStaticsProblem3D(mesh *M, triplet *stiff, csr *K_csr, double **b, double **dirichletValues, double **x, double F_magnitude, double azimuth_deg, double elevation_deg, int gravity);
void ExecuteStaticsAnalysis3D(mesh *M, csr *K_csr, double *b, double *x, int nDOF, int solverType, double tolerance);
void CleanupStaticsProblem3D(triplet *stiff, csr *K_csr, double *b, double *dirichletValues, double *x);
/* ------------------------------------------------------------------------------------------------- */
void InitializeStaticsProblem3D_Orthotropic(mesh *M, triplet *stiff, csr *K_csr, double **b, double **dirichletValues, double **x, double F_magnitude, double azimuth_deg, double elevation_deg, int gravity, double angleX, double angleY, double angleZ);
void ExecuteStaticsAnalysis3D_Orthotropic(mesh *M, csr *K_csr, double *b, double *x, int nDOF, int solverType, double tolerance);
void CleanupStaticsProblem3D_Orthotropic(triplet *stiff, csr *K_csr, double *b, double *dirichletValues, double *x);
/* ------------------------------------------------------------------------------------------------- */
void InitializeDynamicProblem3D(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, int *nDOF, double alpha, double F_magnitude, double azimuth_deg, double elevation_deg, double eps1, double eps2, int gravity, SimulationConfig *config);
void ExecuteDynamicAnalysis3D(mesh *M, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, int nDOF, double omega, double eps1, double eps2, double dt, int nSteps, int Damping, int saveVKT, int saveData, SimulationConfig *config);
void CleanupDynamicProblem3D(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double *u, double *v, double *a, double *F, double *dirichletValues);
/* ------------------------------------------------------------------------------------------------- */
void InitializeDynamicProblem3D_Orthotropic(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, int *nDOF, double alpha, double F_magnitude, double azimuth_deg, double elevation_deg, double eps1, double eps2, int gravity, double angleX, double angleY, double angleZ, int Damping, SimulationConfig *config);
void ExecuteDynamicAnalysis3D_Orthotropic(mesh *M, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, int nDOF, double omega, double eps1, double eps2, double dt, int nSteps, int Damping, int saveVKT, int saveData, SimulationConfig *config);
void CleanupDynamicProblem3D_Orthotropic(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, int Damping, double *u, double *v, double *a, double *F, double *dirichletValues);
/* ------------------------------------------------------------------------------------------------- */
// void InitializeNonlinearElasticityProblem(mesh *M, triplet *stiff, csr *K_csr, double **displacement, double **residual, double **dirichletValues, double F_magnitude, double azimuth_deg, double elevation_deg);
void InitializeNonlinearElasticityProblem(mesh *M, triplet *stiff, csr *K_csr, double **displacement, double **residual, double **f_ext, double **P_int, double **dirichletValues, double F_magnitude, double azimuth_deg, double elevation_deg);
// void ExecuteNonlinearElasticityAnalysis3D(mesh *M, triplet stiff, double *displacement, double *residual, double *dirichletValues, double F_magnitude, double azimuth_deg, double elevation_deg, int maxIter, double tolerance);
void ExecuteNonlinearElasticityAnalysis3D(mesh *M, triplet stiff, double *displacement, double *residual, double *f_ext, double *P_int, double *dirichletValues, double F_magnitude, double azimuth_deg, double elevation_deg, int maxIter, double tolerance);
void CleanupNonlinearElasticityProblem(mesh *M, triplet *stiff, csr *K_csr, double *displacement, double *residual, double *dirichletValues);
/* ------------------------------------------------------------------------------------------------- */
void InitializeNonlinearElasticityProblem_Dynamics(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, double maxIter, double tolerance, double maxNewton, double dt, int nSteps, int useDamping, int saveVibration, int saveData, int *nDOF, double alpha, double F_magnitude, double azimuth_deg, double elevation_deg, double eps1, double eps2, int gravity, SimulationConfig *config);
void CleanupNonlinearElasticityProblem_Dynamics(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double *u, double *v, double *a, double *F, double *dirichletValues);
/* ------------------------------------------------------------------------------------------------- */

// External dependencies
double f(double x, double y);
double f3D(double x, double y, double z);

double psi(double x, double y, double nx, double ny);
double psi3D(double x, double y, double z, double nx, double ny, double nz);

double ur(double x, double y, double nx, double ny);
double ur3D(double x, double y, double z, double nx, double ny, double nz);

double psi_force(double x, double y, double nx, double ny, double angle_degrees, double *force_x, double *force_y);
double psi_force_3D(double x, double y, double z, double nx, double ny, double nz, double F_magnitude, double azimuth_deg, double elevation_deg, double *force_x, double *force_y, double *force_z);

void Mesh_Discretize_Poisson(mesh *this, double f(double, double), triplet *p, double *b, double *dirichletValues, double (*psi)(double, double, double, double), double (*ur)(double, double, double, double), IntegrationType integrationType);
void Mesh_Discretize_Poisson_3D(mesh *this, double (*f)(double, double, double), triplet *p, double *b, double *dirichletValues, double (*psi)(double, double, double, double, double, double), double (*ur)(double, double, double, double, double, double));

void Mesh_Discretize_Statics(mesh *this, double *b, triplet *p, double *dirichletValues, double (*psi_force)(double, double, double, double, double, double *, double *), double angle_degrees, int gravity);
void Mesh_Discretize_Statics3D(mesh *this, double *b, triplet *p, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, int gravity, double angleX, double angleY, double angleZ);
void Mesh_Discretize_Statics3D_Orthotropic(mesh *this, double *b, triplet *p, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, int gravity, double angleX, double angleY, double angleZ);

void Mesh_Discretize_Dynamics(mesh *this, double *b, triplet *mass, triplet *p, triplet *dumping, double *dirichletValues, double (*psi_force)(double, double, double, double, double, double *, double *), double angle_degrees, const double eps1, const double eps2, int gravity);
void Mesh_Discretize_Dynamics3D(mesh *this, double *b, triplet *mass, triplet *p, triplet *damping, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, const double eps1, const double eps2, int gravity, double angleX, double angleY, double angleZ);
void Mesh_Discretize_Dynamics3D_Orthotropic(mesh *this, double *b, triplet *mass, triplet *p, triplet *damping, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, const double eps1, const double eps2, int gravity, double angleX, double angleY, double angleZ);

void Mesh_Discretize_Mass_Matrix(mesh *this, triplet *mass);

void NewmarkIntegrate(mesh *this, csr *M_csr, csr *K_csr, csr *D_csr , double *F, double *u, double *v, double *a, double dt, int nSteps, int useDamping, int saveVibration, int saveData, SimulationConfig *config);
void NewmarkIntegrate_Nonlinear(mesh *M, csr *M_csr, csr *D_csr, double *F_ext, double *u, double *v, double *a, double dt, int nSteps, int maxNewton, double tol, int useDamping, int saveVibration, int saveData, SimulationConfig *config, double (*psi_force_3D)(double,double,double,double,double,double,double,double,double,double*,double*,double*), double *dirichletValues, double F_mag, double azimuth_deg, double elevation_deg);
#endif
