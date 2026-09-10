#pragma once

#include "Triplet_1.h"
#include "vector_1.h"

#define Max(a, b) (((a) > (b)) ? (a) : (b))

typedef struct {
    int     *PI;    // ukazatel na radky
    int     *J;     // indexy sloupcu
    double  *VAL;   // prvky matice
    int     n;      // rozmer matice
    int     nz;     // pocet nenulovych prvku
    int     nalloc; // alokovana pamet
} sparse;

typedef struct {
    int     *PI;   
    int     *J;     
    double  *VAL;   // realna cast
    double  *IVAL;  // imaginarni cast
    int     n;      
    int     nz;     
    int     nalloc; 
} isparse;

void Sparse_Allocate(sparse* this, int n, int nalloc);
void Sparse_Free(sparse* this);
void Sparse_Create(sparse* this, triplet_1* tri);
void Sparse_Multiply(sparse* A, vector x, vector Ax);
void Sparse_Multiply_Konst(sparse* A, double konst);
void Sparse_Add(sparse* A, double konst, sparse* M, sparse* AM, vector isfixed);
void Sparse_Add2(sparse* A, double konst, sparse* B, sparse* AB);
void Sparse_Copy(sparse* A, sparse* B);
void Sparse_GaussSeidel(sparse* A, vector b, vector x, int niter, double err);
void Sparse_SOR(sparse* A, vector b, vector x, int niter, double err, double omega);
void Sparse_SteepestDescentMethod(sparse* A, vector b, vector x, int niter, double err);
void Sparse_ConjugateGradientMethod(sparse* A, vector b, vector x, int niter, double err);
void Sparse_PreconditionedConjugateGradientMethod(sparse* A, vector b, vector x, int niter, double err);
void Sparse_Rezidiuum(sparse* A, vector b, vector x, vector rez);
double Sparse_GetRowRezidiuum(sparse* A, vector b, vector x, int i, double* aii);
void Sparse_SolveUMFPACK(sparse *A, vector b, vector x);
void Sparse_SolveUMFPACK_IMAG(sparse *A, vector A_IM, vector b, vector b_IM, vector x, vector x_IM);
void Sparse_Save(const sparse* this, const char* fname);

void Sparse_Imag_Allocate(isparse* this, int n, int nalloc);
void Sparse_Imag_Free(isparse* this);
// void Sparse_Imag_Create(isparse* this, itriplet* tri);
void Sparse_Imag_Multiply_Konst(isparse* A, double konst);
void Sparse_Imag_Add(isparse* A, double konst, isparse* M, isparse* AM, vector isfixed);
void Sparse_Imag_Copy(isparse* A, isparse* B);
void Sparse_Imag_SolveUMFPACK(isparse *A, vector b, vector ib, vector x, vector ix);
void Sparse_Imag_Save(const isparse* this, const char* fname);