#ifndef NONLINEAR_H
#define NONLINEAR_H

#include "csr.h"
#include "mesh.h"
#include "triplet.h"

// void solveNonLinear(mesh *M, double *u_global, triplet *stiff, double *residual, double *dirichletValues, int maxIter, double tolerance, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg);
void solveNonLinear(mesh *M, double *u_global, triplet *stiff, double *residual, double *f_ext, double *P_int, double *dirichletValues, int maxIter, double tolerance, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg);
void runNonlinearWithLoadSteps(mesh *M, double *u_global, triplet *stiff, double *residual, double *f_ext, double *P_int, double *dirichletValues, int maxIter, double tolerance, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_total, double azimuth_deg, double elev_deg, int nSteps);

// External dependencies
void initializeMaterialProperties(mesh *M, MaterialProperties *materials);
void computeResidualForceVector(mesh *this, double *u_global, double *residual, MaterialProperties *materials);
void computeTangentStiffnessMatrix(mesh *this, int elementIndex, double localK[24][24], MaterialProperties *materials, double *u_global);

void computeExternalForces(mesh *M, double *f_ext, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg);
void ApplyDirichletBCsFromFaces3D(mesh *M, triplet *stiff, double *residual, double *dirichletValues);

// void Mesh_Discretize_NonlinearStatics3D(mesh *M, triplet *stiff, double *residual, double *u_global, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg);
void Mesh_Discretize_NonlinearStatics3D(mesh *M, triplet *stiff, double *residual, double *f_ext, double *P_int, double *u_global, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg);
#endif // NONLINEAR_H
