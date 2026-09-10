#ifndef SAVING_H
#define SAVING_H

#include <stdio.h>
#include "csr.h"  // For csr matrix definition
#include "mesh.h" // For mesh definition

void WriteTripleteMatrixToFile(const triplet *p, const char *filename);
void saveSparseMatrixToFile(const csr *A_csr, const char *filename);

void exportNumericalSolutionToVTK(const mesh *M, const double *numericalSolution, const char *filename);
void exportNumericalSolutionToVTK_3D(const mesh *M, const double *numericalSolution, const char *filename);

void exportDisplacementToVTK(const mesh *M, const double *displacement, const char *filename);
void exportDisplacementToVTK_3D(const mesh *M, const double *displacement, const char *filename);

void saveVectorToFile(double *vector, int size, int dimensions, const char *filename);

void saveTimeStepData(const char *filename, double *displacement, double *velocity, int NPoints, double time);
void saveTimeStepData3D(const char *filename, double *displacement, double *velocity, int NPoints, double time);

void saveTimeSeriesData(FILE *time_series_file, double time, double *u, double *v, double *a, int nDOF);
void saveTimeSeriesData3D(FILE *time_series_file, double time, double *u, double *v, double *a, int nDOF);

void saveIterationValue(const char *filename, int iteration, double value);

#endif // SAVING_H
