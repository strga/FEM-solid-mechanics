#ifndef POISSON_H
#define POISSON_H

#include "mesh.h"

void CalculateAnalyticalGradient(double x, double y, double *grad_x, double *grad_y);
void CalculateAnalyticalSolution(mesh *M, double *analyticalSolution);
double AnalyticalSolution(double x, double y); 
void Poisson_ComputeNorms(mesh *M, double *numericalSolution);

void CalculateAnalyticalGradient3D(double x, double y, double z, double *grad_x, double *grad_y, double *grad_z);
void CalculateAnalyticalSolution3D(mesh *M, double *analyticalSolution3D);
double AnalyticalSolution3D(double x, double y, double z);
void Poisson_ComputeNorms3D(mesh *M, double *numericalSolution3D);

#endif // POISSON_H
