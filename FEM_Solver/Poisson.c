#include <math.h>
#include <stdio.h>
#include "Poisson.h"
#include "mesh.h"
#include "constants.h"
// Analytical gradient computation for Poisson
void CalculateAnalyticalGradient(double x, double y, double *grad_x, double *grad_y) {
    *grad_x = M_PI * cos(M_PI * x) * sin(M_PI * y);
    *grad_y = M_PI * sin(M_PI * x) * cos(M_PI * y);
}

double AnalyticalSolution(double x, double y) {
    return sin(M_PI * x) * sin(M_PI * y);
}

// Analytical solution for Poisson
void CalculateAnalyticalSolution(mesh *M, double *analyticalSolution) {
    for (int i = 0; i < M->NPoints; ++i) {
        double x = M->x[i];
        double y = M->y[i];
        analyticalSolution[i] = AnalyticalSolution(x, y);
    }
}

// Computes and prints norms for Poisson problem
void Poisson_ComputeNorms(mesh *M, double *numericalSolution) {

    double L2NormError = computeL2Error(M, numericalSolution, AnalyticalSolution);
    double H1SeminormError = computeH1SeminormError(M, numericalSolution, AnalyticalSolution, CalculateAnalyticalGradient);

    printf("NELEMENTS:%d \t Points:%d\n", M->NElements, M->NPoints);
    printf("Step size: h = 1 / %d\n", (int)(sqrt(M->NPoints + 1.0) - 1.0));
    printf("L2 Norm of the solution error: %g\n", L2NormError);
    printf("H1 Seminorm of the gradient error: %g\n", H1SeminormError);
}

// ------------------------------------------------------------------------------------- //
// 3D version
// ------------------------------------------------------------------------------------- //

// Analytical gradient computation for Poisson
void CalculateAnalyticalGradient3D(double x, double y, double z, double *grad_x, double *grad_y, double *grad_z) {
    *grad_x = M_PI * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
    *grad_y = M_PI * sin(M_PI * x) * cos(M_PI * y) * sin(M_PI * z);
    *grad_z = M_PI * sin(M_PI * x) * cos(M_PI * y) * cos(M_PI * z);
}

double AnalyticalSolution3D(double x, double y, double z) {
    return sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
}

// Analytical solution for Poisson
void CalculateAnalyticalSolution3D(mesh *M, double *analyticalSolution3D) {
    for (int i = 0; i < M->NPoints; ++i) {
        double x = M->x[i];
        double y = M->y[i];
        double z = M->z[i];
        analyticalSolution3D[i] = AnalyticalSolution3D(x, y, z);
    }
}

// Computes and prints norms for Poisson problem
void Poisson_ComputeNorms3D(mesh *M, double *numericalSolution3D) {

    double L2NormError3D = computeL2Error3D(M, numericalSolution3D, AnalyticalSolution3D);
    double H1SeminormError3D = computeH1SeminormError3D(M, numericalSolution3D, AnalyticalSolution3D, CalculateAnalyticalGradient3D);

    printf("NELEMENTS:%d \t Points:%d\n", M->NElements, M->NPoints);
    printf("Step size (here with cube root): h = 1 / %d\n", (int)(cbrt(M->NPoints + 1.0) - 1.0));
    printf("L2 Norm of the solution error: %g\n", L2NormError3D);
    printf("H1 Seminorm of the gradient error: %g\n", H1SeminormError3D);
}

