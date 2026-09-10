//#define DEBUG
#include <umfpack.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sparse_1.h"
#include "Triplet_1.h"


void Sparse_Allocate(sparse* this, int n, int nalloc) {
    this->PI = malloc((n + 1) * sizeof(int));
    if (this->PI == NULL) {
        printf("\nError[Sparse_Allocate]: Nedostatek pameti pro sparse.PI\n");
        exit(1);
    }

    this->J = malloc(nalloc * sizeof(int));
    if (this->J == NULL) {
        printf("\nError[Sparse_Allocate]: Nedostatek pameti pro sparse.J\n");
        exit(1);
    }

    this->VAL = malloc(nalloc * sizeof(double));
    if (this->VAL == NULL) {
        printf("\nError[Sparse_Allocate]: Nedostatek pameti pro sparse.VAL\n");
        exit(1);
    }

    this->n      = n;
    this->nz     = 0;
    this->nalloc = nalloc;
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Free(sparse* this) {
    if (this->nalloc > 0) {
        free(this->PI);
        free(this->J);
        free(this->VAL);

        this->nalloc = 0;
        this->nz     = 0;
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Create(sparse* this, triplet_1* tri) {
    Triplet_1_Quicksort(tri, 0, tri->nz);
    Triplet_1_Unique(tri);
    Sparse_Allocate(this, tri->n, tri->nz);

    int i, j = 0;
#define DEBUG1
#ifdef DEBUG1
    printf("SPARSE_CREATE:\n");
    printf("tri.n = %d\n", tri->n);
    printf("tri.nz = %d\n", tri->nz);
    printf("sparse.n = %d\n", this->n);
    printf("\n");
    getchar();
    // printf("\ntriplet:\n");
    // for (i = 0; i < tri->nz; i++) {
    //     printf("%d %d %lf\n", tri->I_re[i], tri->J[i], tri->VAL[i]);
    // }
#endif // DEBUG

    this->PI[0] = 0;

    for (i = 0; i < this->n; i++) {
        while (tri->I_re[j] == i && j < tri->nz) {
            this->J[j]   = tri->J[j];
            this->VAL[j] = tri->VAL[j];
            j++;
        }
        this->PI[i + 1] = j;
    }
    this->nz = j;

#ifdef DEBUG
    printf("\nsparse.nz = %d\n", this->nz);

    printf("\npole PI pro format sparse:\n");
    for (i = 0; i < this->n + 1; i++) {
        printf("PI[%d] = %d \n", i, this->PI[i]);
    }

    printf("\nsparse matice:\n");
    for (i = 0; i < this->n; i++) {
        for (j = this->PI[i]; j < this->PI[i + 1]; j++) {
            printf("%d %d %lf\n", i, this->J[j], this->VAL[j]);
        }
    }

    printf("\nreziduum:\n");
#endif // DEBUG
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Multiply_Konst(sparse* A, double konst) {
    int i, j;
    
    for (i = 0; i < A->n; i++) {
        for (j = A->PI[i]; j < A->PI[i+1]; j++) {
            A->VAL[j] *= konst;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Multiply(sparse* A, vector x, vector Ax) {
    int k, i, j;
    double aij;

    for (i = 0; i < A->n; i++) {
        Ax[i] = 0;
        for (k = A->PI[i]; k < A->PI[i + 1]; k++) {
            aij    = A->VAL[k];
            j      = A->J[k];
            Ax[i] += aij * x[j];
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Add(sparse* A, double konst, sparse* M, sparse* AM, vector isfixed) {
    int i, j, k, m1, m2;
    m1 = m2 = 0;
    
    if (A->n != M->n) {
        printf("\nError[Sparse_Add]: Rozmery matic se neshoduji\n");
        printf("A.n = %d a M.n = %d\n", A->n, M->n);
        exit(1);
    }

    AM->PI[0] = 0;
    
    for (i = 0; i < A->n; i++) {
        if (isfixed[i] == 1) {
            AM->J[m1]    = A->J[A->PI[i]];
            AM->VAL[m1]  = A->VAL[A->PI[i]];
            AM->PI[m2+1] = AM->PI[m2]+1;
            m1++; m2++;
        }

        else {
            k = M->PI[i];

            for (j = A->PI[i]; j < A->PI[i+1]; j++) {
                while (A->J[j] >= M->J[k] && k < M->PI[i+1]) {

                    if (A->J[j] == M->J[k]) {
                        AM->J[m1]   = A->J[j];
                        AM->VAL[m1] = A->VAL[j] + konst * M->VAL[k];
                        m1++;
                    }

                    else {
                        if(isfixed[M->J[k]] != 1) {
                            AM->J[m1]   = M->J[k];
                            AM->VAL[m1] = konst * M->VAL[k];
                            m1++;
                        }
                    }
                    k++;
                }
            }
            m2++;
            AM->PI[m2] = m1;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Add2(sparse* M, double konst, sparse* D, sparse* MD) {
    int i, j, k, MDpos, MDnz;
    MDpos = MDnz = 0;

    if (M->n != D->n) {
        printf("\nError[Sparse_Add2]: Rozmery matic se neshoduji\n");
        exit(1);
    }

    MD->PI[0] = 0;

    for (i = 0; i < M->n; i++) {
        k = D->PI[i];
        
        for (j = M->PI[i]; j < M->PI[i + 1]; j++) {
            
            if (k == D->PI[i+1]) {
                MD->J[MDnz]   = M->J[j];
                MD->VAL[MDnz] = M->VAL[j];
                MDnz++;
            }

            if (M->J[j] < D->J[k] && k < D->PI[i+1]) {
                MD->J[MDnz]   = M->J[j];
                MD->VAL[MDnz] = M->VAL[j];
                MDnz++;
            }

            while (M->J[j] >= D->J[k] && k < D->PI[i+1]) {
                if (M->J[j] == D->J[k]) {
                    MD->J[MDnz]   = M->J[j];
                    MD->VAL[MDnz] = M->VAL[j] + konst * D->VAL[k];
                    MDnz++; 
                }
                else{
                    MD->J[MDnz]   = D->J[k];
                    MD->VAL[MDnz] = konst * D->VAL[k];
                    MDnz++; 
                }
                k++;
            }

            

            
                

        }
        MDpos++;
        MD->PI[MDpos] = MDnz;
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Copy(sparse* A, sparse* B) {
    B->PI[0] = 0;

    for (int i = 0; i < A->n; i++) {
        for (int j = A->PI[i]; j < A->PI[i + 1]; j++) {
            B->J[j]   = A->J[j];
            B->VAL[j] = A->VAL[j];
        }
        B->PI[i+1] = A->PI[i+1];
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_GaussSeidel(sparse* A, vector b, vector x, int niter, double err) {
    int iter, i;
    double aii, rez, rezi;

    for (iter = 0; iter < niter; iter++) {
        rez = 0;
        for (i = 0; i < A->n; i++) {
            rezi  = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
            x[i] += rezi / aii;
            rez   = Max(rez, fabs(rezi));
        }

#ifdef DEBUG
        printf("rez[%d] = %e\n", iter, rez);
#endif // DEBUG

        if (rez < err) {
            break;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_SOR(sparse * A, vector b, vector x, int niter, double err, double omega) {
    int iter, i;
    double aii, rez, rezi;

    for (iter = 0; iter < niter; iter++) {
        rez = 0;
        for (i = 0; i < A->n; i++) {
            rezi  = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
            x[i] += omega * rezi / aii;
            rez   = Max(rez, fabs(rezi));
        }

#ifdef DEBUG
        printf("rez[%d] = %e\n", iter, rez);
#endif // DEBUG

        if (rez < err) {
            break;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_SteepestDescentMethod(sparse* A, vector b, vector x, int niter, double err) {
    int iter, i;
    double aii, rez, alfa, pr, pAp;
    vector p, Ap;

    p  = Vector_Allocate(A->n);
    Ap = Vector_Allocate(A->n);

    for (i = 0; i < A->n; i++) {
        p[i] = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
    }

    for (iter = 0; iter < niter; iter++) {
        rez = pr = pAp = 0;
        Sparse_Multiply(A, p, Ap);

        for (i = 0; i < A->n; i++) {
            pr  += p[i] * p[i];
            pAp += p[i] * Ap[i];
        }
        alfa = pr / pAp;

        for (i = 0; i < A->n; i++) {
            x[i] += alfa * p[i];
            p[i] -= alfa * Ap[i];
            rez   = Max(rez, fabs(p[i]));
        }

#ifdef DEBUG
        printf("rez[%d] = %e\n", iter, rez);
#endif // DEBUG

        if (rez < err) {
            break;
        }
    }

    Vector_Free(p);
    Vector_Free(Ap);
}
/* ----------------------------------------------------------------------------------- */
void Sparse_ConjugateGradientMethod(sparse * A, vector b, vector x, int niter, double err) {
    int iter, i;
    double aii, rez, alfa, pr, pAp;
    vector p, Ap;
    double beta, rAp, rr, rrnew;
    vector rezi;

    p    = Vector_Allocate(A->n);
    Ap   = Vector_Allocate(A->n);
    rezi = Vector_Allocate(A->n);

    for (i = 0; i < A->n; i++) {
        p[i] = rezi[i] = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
    }

    for (iter = 0; iter < niter; iter++) {
        rez = pr = pAp = rr = rrnew = rAp = 0;

        Sparse_Multiply(A, p, Ap);
        for (i = 0; i < A->n; i++) {
            pr  += p[i] * rezi[i];  
            pAp += p[i] * Ap[i];
        }
        alfa = pr / pAp;     

        for (i = 0; i < A->n; i++) {
            x[i]    += alfa * p[i];
            rezi[i] -= alfa * Ap[i];
            rez      = Max(rez, fabs(rezi[i]));
        }

        for (i = 0; i < A->n; i++) {
            rAp += rezi[i] * Ap[i];
        }
        beta = -rAp / pAp;   

        for (i = 0; i < A->n; i++) {
            p[i] = rezi[i] + beta * p[i];
        }

#ifdef DEBUG
        printf("rez[%d] = %e\n", iter, rez);
#endif // DEBUG

        if (rez < err) {
            break;
        }
    }

    Vector_Free(p);
    Vector_Free(Ap);
    Vector_Free(rezi);
}
/* ----------------------------------------------------------------------------------- */
void Sparse_PreconditionedConjugateGradientMethod(sparse* A, vector b, vector x, int niter, double err) {
    int iter, i, j;
    double aii, rez, alfa, pAp;
    vector p, Ap;
    double beta;
    vector rezi;
    double rz, rznew;
    vector M, z;

    p    = Vector_Allocate(A->n);
    Ap   = Vector_Allocate(A->n);
    rezi = Vector_Allocate(A->n);
    M    = Vector_Allocate(A->n);
    z    = Vector_Allocate(A->n);

    for (i = 0; i < A->n; i++) {
        for (j = A->PI[i]; j < A->PI[i + 1]; j++) {
            if (i == A->J[j]) {
                M[i] = A->VAL[j];
            }
        }
    }

    for (i = 0; i < A->n; i++) {
        rezi[i] = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
        p[i] = z[i] = rezi[i] / M[i];
    }

    for (iter = 0; iter < niter; iter++) {
        rez = pAp = rz = rznew = 0;

        Sparse_Multiply(A, p, Ap);
        for (i = 0; i < A->n; i++) {
            rz  += rezi[i] * z[i];
            pAp += p[i] * Ap[i];
        }
        alfa = rz / pAp;

        for (i = 0; i < A->n; i++) {
            x[i]    += alfa * p[i];
            rezi[i] -= alfa * Ap[i];
            rez      = Max(rez, fabs(rezi[i]));
            z[i]     = rezi[i] / M[i];
        }

        for (i = 0; i < A->n; i++) {
            rznew += rezi[i] * z[i];
        }
        beta = rznew / rz;

        for (i = 0; i < A->n; i++) {
            p[i] = z[i] + beta * p[i];
        }

#ifdef DEBUG
        printf("rez[%d] = %e\n", iter, rez);
#endif // DEBUG
//printf("rez[%d] = %e\n", iter, rez);

        if (rez < err) {
            break;
        }
    }

    Vector_Free(p);
    Vector_Free(Ap);
    Vector_Free(rezi);
    Vector_Free(M);
    Vector_Free(z);
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Rezidiuum(sparse* A, vector b, vector x, vector rez) {
    int i;
    double aii;

    for (i = 0; i < A->n; i++)
        rez[i] = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
}
/* ----------------------------------------------------------------------------------- */
double Sparse_GetRowRezidiuum(sparse* A, vector b, vector x, int i, double* aii) {
    int k, j;
    double rezi, aij;

    rezi = b[i];

    for (k = A->PI[i]; k < A->PI[i + 1]; k++) {
        aij   = A->VAL[k];
        j     = A->J[k];
        rezi -= aij * x[j];

        if (i == j)
            *aii = aij;
    }
    return rezi;
}
/* ----------------------------------------------------------------------------------- */
void Sparse_SolveUMFPACK(sparse *A, vector b, vector x) {
    void   *SolverSymbolic;
    void   *SolverNumeric;
    void   *SolverControl = NULL;
    void   *SolverInfo    = NULL;

    int    res1, res2, res3;

    res1 = umfpack_di_symbolic(A->n, A->n, A->PI, A->J, A->VAL, &SolverSymbolic, SolverControl, SolverInfo);
    res2 = umfpack_di_numeric(A->PI, A->J, A->VAL, SolverSymbolic, &SolverNumeric, SolverControl, SolverInfo);

    res3 = umfpack_di_solve (UMFPACK_At, A->PI, A->J, A->VAL, x, b, SolverNumeric, SolverControl, SolverInfo);

    umfpack_di_free_symbolic(&SolverSymbolic);
    umfpack_di_free_numeric(&SolverNumeric);

    if (res1 != 0 || res2 != 0 || res3 != 0) {
        printf("\nError[Sparse_SolveUMFPACK]: res1 = %d, res2 = %d, res3 = %d\n", res1, res2, res3);
        exit(1);
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_SolveUMFPACK_IMAG(sparse *A, vector A_IM, vector b, vector b_IM, vector x, vector x_IM) {
    void   *SolverSymbolic;
    void   *SolverNumeric;
    void   *SolverControl = NULL;
    void   *SolverInfo    = NULL;

    int    res1, res2, res3;

    res1 = umfpack_zi_symbolic (A->n, A->n, A->PI, A->J, A->VAL, A_IM, &SolverSymbolic, SolverControl, SolverInfo);
    res2 = umfpack_zi_numeric (A->PI, A->J, A->VAL, A_IM, SolverSymbolic, &SolverNumeric, SolverControl, SolverInfo);

    res3 = umfpack_zi_solve (UMFPACK_Aat, A->PI, A->J, A->VAL, A_IM, x, x_IM, b, b_IM, SolverNumeric, SolverControl, SolverInfo);

    umfpack_zi_free_symbolic (&SolverSymbolic);
    umfpack_zi_free_numeric (&SolverNumeric);

if (res1 != 0 || res2 != 0 || res3 != 0) {
        printf("\nError[Sparse_SolveUMFPACK_IMAG]: res1 = %d, res2 = %d, res3 = %d\n", res1, res2, res3);
        exit(1);
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Save(const sparse* this, const char* fname) {
    int   i, j;
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Sparse_Save]: Nepodarilo se otevrit soubor pro SPARSE matici\n");
        exit(1);
    }

    for (i = 0; i < this->n; i++) {
        for (j = this->PI[i]; j < this->PI[i + 1]; j++) {
            fprintf(fid, "%d %d %e\n", i + 1, this->J[j] + 1, this->VAL[j]);
        }
    }

    fclose(fid);
}


void Sparse_Imag_Allocate(isparse* this, int n, int nalloc) {
    this->PI = malloc((n + 1) * sizeof(int));
    if (this->PI == NULL) {
        printf("\nError[Sparse_Imag_Allocate]: Nedostatek pameti pro sparse.PI\n");
        exit(1);
    }

    this->J = malloc(nalloc * sizeof(int));
    if (this->J == NULL) {
        printf("\nError[Sparse_Imag_Allocate]: Nedostatek pameti pro sparse.J\n");
        exit(1);
    }

    this->VAL = malloc(nalloc * sizeof(double));
    if (this->VAL == NULL) {
        printf("\nError[Sparse_Imag_Allocate]: Nedostatek pameti pro sparse.VAL\n");
        exit(1);
    }

    this->IVAL = malloc(nalloc * sizeof(double));
    if (this->IVAL == NULL) {
        printf("\nError[Sparse_Imag_Allocate]: Nedostatek pameti pro sparse.IVAL\n");
        exit(1);
    }

    this->n      = n;
    this->nz     = 0;
    this->nalloc = nalloc;
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_Free(isparse* this) {
    if (this->nalloc > 0) {
        free(this->PI);
        free(this->J);
        free(this->VAL);
        free(this->IVAL);

        this->nalloc = 0;
        this->nz     = 0;
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_Create(isparse* this, itriplet* tri) {
    Triplet_Imag_Quicksort(tri, 0, tri->nz);
    Triplet_Imag_Unique(tri);
    Sparse_Imag_Allocate(this, tri->n, tri->nz);

    int i, j = 0;

#ifdef DEBUG
    printf("\ntri.n = %d\n", tri->n);
    printf("tri.nz = %d\n", tri->nz);
    printf("sparse.n = %d\n", this->n);
    printf("\nitriplet:\n");
    for (i = 0; i < tri->nz; i++) {
        printf("%d %d %lf %lf\n", tri->I_re[i], tri->J[i], tri->VAL[i], tri->IVAL[i]);
    }
#endif // DEBUG

    this->PI[0] = 0;

    for (i = 0; i < this->n; i++) {
        while (tri->I_im[j] == i + 1 && j < tri->nz) {
            this->J[j]   = tri->J[j] - 1;
            this->VAL[j] = tri->VAL[j];
            this->IVAL[j] = tri->IVAL[j];
            j++;
        }
        this->PI[i + 1] = j;
    }
    this->nz = j;

#ifdef DEBUG
    printf("\nsparse.nz = %d\n", this->nz);

    printf("\npole PI pro format sparse:\n");
    for (i = 0; i < this->n + 1; i++) {
        printf("PI[%d] = %d \n", i, this->PI[i]);
    }

    printf("\nisparse matice:\n");
    for (i = 0; i < this->n; i++) {
        for (j = this->PI[i]; j < this->PI[i + 1]; j++) {
            printf("%d %d %lf %lf\n", i, this->J[j], this->VAL[j], this->IVAL[j]);
        }
    }

    printf("\nreziduum:\n");
#endif // DEBUG
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_Multiply_Konst(isparse* A, double konst) {
    int i, j;
    
    for (i = 0; i < A->n; i++) {
        for (j = A->PI[i]; j < A->PI[i+1]; j++) {
            A->IVAL[j] *= konst;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_Add(isparse* A, double konst, isparse* M, isparse* AM, vector isfixed) {
    int i, j, k, m1, m2;
    m1 = m2 = 0;
    
    if (A->n != M->n) {
        printf("\nError[Sparse_Imag_Add]: Rozmery matic se neshoduji\n");
        exit(1);
    }

    AM->PI[0] = 0;
    
    for (i = 0; i < A->n; i++) {
        if (isfixed[i] == 1) {
            AM->J[m1]    = A->J[A->PI[i]];
            AM->VAL[m1]  = A->VAL[A->PI[i]];
            AM->IVAL[m1] = A->IVAL[A->PI[i]];
            AM->PI[m2+1] = AM->PI[m2] + 1;
            m1++; m2++;
        }

        else {
            k = M->PI[i];

            for (j = A->PI[i]; j < A->PI[i+1]; j++) {
                while (A->J[j] >= M->J[k] && k < M->PI[i+1]) {
                    
                    if (A->J[j] == M->J[k]) {
                        AM->J[m1]    = A->J[j];
                        AM->VAL[m1]  = A->VAL[j] + konst * M->VAL[k];
                        AM->IVAL[m1] = A->IVAL[j];
                        m1++;
                    }

                    k++;
                }
            }

            m2++;
            AM->PI[m2] = m1;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_Copy(isparse* A, isparse* B) {
    B->PI[0] = 0;

    for (int i = 0; i < A->n; i++) {
        for (int j = A->PI[i]; j < A->PI[i + 1]; j++) {
            B->J[j]    = A->J[j];
            B->VAL[j]  = A->VAL[j];
            B->IVAL[j] = A->IVAL[j];
        }
        B->PI[i+1] = A->PI[i+1];
    }
/*
    printf("\npole PI pro format sparse:\n");
    for (int i = 0; i < B->n + 1; i++) {
        printf("B.PI[%d] = %d \t A.PI[%d] = %d \n", i, B->PI[i], i, A->PI[i]);
    }
*/
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_SolveUMFPACK(isparse *A, vector b, vector ib, vector x, vector ix) {
    void   *SolverSymbolic;
    void   *SolverNumeric;
    void   *SolverControl = NULL;
    void   *SolverInfo    = NULL;

    int    res1, res2, res3;

    res1 = umfpack_zi_symbolic (A->n, A->n, A->PI, A->J, A->VAL, A->IVAL, &SolverSymbolic, SolverControl, SolverInfo);
    res2 = umfpack_zi_numeric (A->PI, A->J, A->VAL, A->IVAL, SolverSymbolic, &SolverNumeric, SolverControl, SolverInfo);

    res3 = umfpack_zi_solve (UMFPACK_Aat, A->PI, A->J, A->VAL, A->IVAL, x, ix, b, ib, SolverNumeric, SolverControl, SolverInfo);

    umfpack_zi_free_symbolic (&SolverSymbolic);
    umfpack_zi_free_numeric (&SolverNumeric);

if (res1 != 0 || res2 != 0 || res3 != 0) {
        printf("\nError[Sparse_SolveUMFPACK_IMAG]: res1 = %d, res2 = %d, res3 = %d\n", res1, res2, res3);
        exit(1);
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Imag_Save(const isparse* this, const char* fname) {
    int   i, j;
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Sparse_Imag_Save]: Nepodarilo se otevrit soubor pro SPARSE matici\n");
        exit(1);
    }

    for (i = 0; i < this->n; i++) {
        for (j = this->PI[i]; j < this->PI[i + 1]; j++) {
            fprintf(fid, "%d\t%d\t\t%e\t%e\n", i + 1, this->J[j] + 1, this->VAL[j], this->IVAL[j]);
        }
    }

    fclose(fid);
}
