#ifndef CSR_H
#define CSR_H

#include "triplet.h"

//#define max(a, b) (((a) > (b)) ? (a) : (b))

// csr - compressed sparse row (Yale)
typedef struct{

    int *PI;
    int *J;

    double *VAL;

    int n;
    int nz;
    int nalloc;
}csr;

void CSR_Allocate( csr* this, int n, int nalloc);
void CSR_Free( csr* this); //int n, int nalloc
void PrintCSR(const char *label, csr *this);

void CSR_InitFromTriplet( csr* this, triplet* tri);

void CSR_Multiply(csr *A, double *x, double *Ax);
void CSR_AddScaled(csr *A, csr *B, csr *result, double scalarA, double scalarB);
void CSR_Add(csr *A, csr *B, csr *result);
void CSR_ScalarMultiply(csr *A, double scalar, csr *result);

void CSR_Consolidate(csr *A);
double CSR_RowResiduum(csr* A, const double *b, const double *x, int i, double *aii);

void Gauss_Seidel(csr *A, const double *b, double *x, int max_iter, double tolerance);
void SOR(csr *A, const double *b, double *x, int max_iter, double tolerance, double omega);
void ConjugateGradient(csr* A, const double* b, double* x, int max_iter, double tolerance);
void CSR_Solve(csr *A, const double *b, double *x);

void scaleMatricesForDirichlet(mesh *this, csr *M, csr *D, double omega, double eps1, double eps2, int nDOF_per_node, int useDamping);
#endif