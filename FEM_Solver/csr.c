#include "csr.h"
#include "mesh.h"
#include "constants.h"
#include "gelement3D.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <umfpack.h>
// #define PRNT_STEPS
// #define TIMER

// #define PRINT_INIT_TRIPLET
// #define DEBUG_print_CSR_TRIPLET
// #define DEBUG_CSR_MULTIPLY
// #define DEBUG_print_SOR
// #define DEBUG_print_RowResiduum
// #define DEBUG_print_AddScaled

// #define DEBUG_print_all
// #define quit

void CSR_Allocate(csr* this, int n, int nalloc)
{
    #ifdef DEBUG_print_all
    printf("Entering CSR_Allocate: n = %d, nalloc = %d\n", n, nalloc);
    #endif

    this->PI = malloc((n + 1) * sizeof(int));
    if (this->PI == NULL)
    {
        printf("Memory Error CSR_Allocate: PI\n");
        exit(EXIT_FAILURE);
    } 
    #ifdef DEBUG_print_all
    else {
        printf("PI allocated successfully.\n");
    }
    #endif

    this->J = malloc(nalloc * sizeof(int));
    if (this->J == NULL)
    {
        printf("Memory Error CSR_Allocate: J\n");
        exit(EXIT_FAILURE);
    }
    #ifdef DEBUG_print_all
    else {
        printf("J allocated successfully.\n");
    }
    #endif

    this->VAL = malloc(nalloc * sizeof(double));
    if (this->VAL == NULL)
    {
        printf("Memory Error CSR_Allocate: VAL\n");
        exit(EXIT_FAILURE);
    } 
    #ifdef DEBUG_print_all
    else {
        printf("VAL allocated successfully.\n");
    }
    #endif

    this->n = n;
    this->nalloc = nalloc;

    #ifdef DEBUG_print_all
    printf("Exiting CSR_Allocate\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_Free(csr* this)
{
    #ifdef DEBUG_print_all
    printf("Entering CSR_Free\n");
    #endif

    if (this != NULL)
    {
        // printf("1\n");
        free(this->PI);
        // printf("2\n");
        free(this->J);
        // printf("3\n");
        free(this->VAL);
        // printf("4\n");

        this->nalloc = 0;
        this->nz = 0;
        this->n = 0;

        #ifdef DEBUG_print_all
        printf("CSR memory freed successfully.\n");
        #endif

    }
    #ifdef DEBUG_print_all
    else {
       
        printf("Error: Attempted to free a null pointer.\n");
    }
    printf("Exiting CSR_Free\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void PrintCSR(const char *label, csr *this) {
    printf("\n%s:\n", label);
    printf("n (number of rows): %d\n", this->n);
    printf("nz (number of non-zero entries): %d\n", this->nz);

    printf("\nRow pointers (PI):\n");
    for (int i = 0; i < this->n + 1; i++) {
        printf("PI[%d] = %d\n", i, this->PI[i]);
    }

    printf("\nColumn indices (J) and Values (VAL):\n");
    for (int i = 0; i < this->n; i++) {
        for (int j = this->PI[i]; j < this->PI[i + 1]; j++) {
            printf("Row %d -> Col %d, VAL = %.10e\n", i, this->J[j], this->VAL[j]);
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_InitFromTriplet(csr* this, triplet* tri)
{
    printf("Triplet_Sort_CSRFromTriplet\n");
    Triplet_Sort1(tri, 0, tri->nz-1); 
    printf("Triplet_Unique_CSRFromTriplet\n");
    Triplet_Unique(tri);
    printf("CSR_Allocate_CSRFromTriplet\n");
    CSR_Allocate(this, tri->n, tri->nz);
    printf("Non-zero(nz) elements: %d \t Total DOFs(n): %d \n", tri->nz, tri->n);

    int i, k = 0;

    #ifdef DEBUG1
    printf("\nCSR_DATA:\n");
    printf("tri.n = %d\n", tri->n);
    printf("tri.nz = %d\n", tri->nz);
    printf("sparse(CSR).n = %d\n", this->n);
    printf("\n");
    // getchar();
    // printf("\ntriplet:\n");
    // for (i = 0; i < tri->nz; i++) {
    //     printf("%d %d %lf\n", tri->tripletI[i], tri->J[i], tri->VAL[i]);
    // }
    #endif

    this->PI[0] = 0;

    for ( i = 0; i < this->n; i++) {
        while (k < tri->nz && tri->tripletI[k] == i) {
            this->J[k] = tri->J[k];
            this->VAL[k] = tri->VAL[k];
            k++;
        }
        this->PI[i + 1] = k;
    }
    this->nz = k;

    #ifdef DEBUG
    printf("\nsparse.nz = %d\n", this->nz);

    printf("\npole PI pro format sparse:\n");
    for (i = 0; i < this->n + 1; i++) {
        printf("PI[%d] = %d \n", i, this->PI[i]);
    }

    printf("\nsparse matice:\n");
    for (i = 0; i < this->n; i++) {
        for (int j = this->PI[i]; j < this->PI[i + 1]; j++) {
            printf("%d %d %lf\n", i, this->J[j], this->VAL[j]);
        }
    }

    printf("\nreziduum:\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_Multiply(csr *A, double *x, double *Ax)
{

    #ifdef DEBUG_CSR_MULTIPLY
    printf("Entering CSR_Multiply\n");
    printf("Matrix dimension: %d\n", A->n);
    printf("Allocated non-zero entries: %d\n", A->nalloc);
    #endif

    int i, j, k;
    double aij;

    for ( i = 0; i < A->n; i++)
    {
        Ax[i] = 0;
        
        #ifdef DEBUG_CSR_MULTIPLY
        printf("Processing row %d\n", i);
        #endif
        
        for ( k = A->PI[i]; k < A->PI[i+1]; k++)
        {
            aij = A->VAL[k];
            j = A->J[k];
            Ax[i] += aij * x[j];

            #ifdef DEBUG_CSR_MULTIPLY
            printf("Contribution to Ax[%d]: %f * x[%d](%f) = %f, Running sum: %f\n", i, aij, j, x[j], aij * x[j], Ax[i]);
            #endif

        }

        #ifdef DEBUG_CSR_MULTIPLY
        printf("Completed row %d, Ax[%d] = %f\n", i, i, Ax[i]);
        #endif
        
    }
    #ifdef DEBUG_CSR_MULTIPLY
    printf("Exiting CSR_Multiply\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_ScalarMultiply(csr *A, double scalar, csr *result) {
    // Allocate result CSR matrix with the same structure as A
    CSR_Allocate(result, A->n, A->nalloc);

    // result->PI = malloc((A->n + 1) * sizeof(int));
    // result->J = malloc(A->nalloc * sizeof(int));
    // result->VAL = malloc(A->nalloc * sizeof(double));

    if (!result->PI || !result->J || !result->VAL) {
        fprintf(stderr, "Memory allocation failed in CSR_ScalarMultiply.\n");
        exit(EXIT_FAILURE);
    }

    result->n = A->n;
    result->nalloc = A->nalloc;
    result->nz = A->nz;

    // Copy row pointers and column indices directly
    memcpy(result->PI, A->PI, (A->n + 1) * sizeof(int));
    memcpy(result->J, A->J, A->nalloc * sizeof(int));

    // Multiply each value in A by scalar and store in result
    for (int i = 0; i < A->nz; i++) {
        result->VAL[i] = scalar * A->VAL[i];
    }

    #ifdef DEBUG_print_all
    printf("CSR_ScalarMultiply: Multiplied matrix by scalar %.5f.\n", scalar);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_Add(csr *A, csr *B, csr *result) {
    assert(A->n == B->n);

    // Allocate result CSR matrix with sufficient space
    CSR_Allocate(result, A->n, A->nalloc + B->nalloc);  // Conservative allocation

    int *tempJ = malloc((A->nz + B->nz) * sizeof(int));
    double *tempVAL = malloc((A->nz + B->nz) * sizeof(double));

    if (!tempJ || !tempVAL) {
        fprintf(stderr, "Memory allocation failed in CSR_Add.\n");
        exit(EXIT_FAILURE);
    }

    int nz = 0;  // Non-zero count
    for (int i = 0; i < A->n; i++) {
        result->PI[i] = nz;

        int a_ptr = A->PI[i];
        int b_ptr = B->PI[i];

        // Merge non-zero elements from both matrices
        while (a_ptr < A->PI[i + 1] || b_ptr < B->PI[i + 1]) {
            if (a_ptr < A->PI[i + 1] && (b_ptr >= B->PI[i + 1] || A->J[a_ptr] < B->J[b_ptr])) {
                // Unique entry from A
                tempJ[nz] = A->J[a_ptr];
                tempVAL[nz] = A->VAL[a_ptr];
                a_ptr++;
            } else if (b_ptr < B->PI[i + 1] && (a_ptr >= A->PI[i + 1] || B->J[b_ptr] < A->J[a_ptr])) {
                // Unique entry from B
                tempJ[nz] = B->J[b_ptr];
                tempVAL[nz] = B->VAL[b_ptr];
                b_ptr++;
            } else {
                // Matching columns, combine values
                tempJ[nz] = A->J[a_ptr];
                tempVAL[nz] = A->VAL[a_ptr] + B->VAL[b_ptr];
                a_ptr++;
                b_ptr++;
            }
            nz++;
        }
    }
    result->PI[A->n] = nz;
    result->nz = nz;

    // Finalize by copying into result->J and result->VAL
    result->J = realloc(tempJ, nz * sizeof(int));
    result->VAL = realloc(tempVAL, nz * sizeof(double));
    if (!result->J || !result->VAL) {
        fprintf(stderr, "Memory reallocation failed in CSR_Add.\n");
        exit(EXIT_FAILURE);
    }

    #ifdef DEBUG_print_all
    printf("CSR_Add: Completed element-wise addition.\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_AddScaled(csr *A, csr *B, csr *result, double scalarA, double scalarB) {
    assert(A->n == B->n);
    int n = A->n;
    result->n = n;
    result->nz = A->nz + B->nz; // adjust as needed
    result->PI = malloc((n + 1) * sizeof(int));
    result->J = malloc(result->nz * sizeof(int));
    result->VAL = malloc(result->nz * sizeof(double));

    #ifdef DEBUG_print_AddScaled
        printf("Starting CSR_AddScaled:\n");
        printf("Matrix dimensions: n = %d\n", n);
        printf("Scalar A: %f, Scalar B: %f\n", scalarA, scalarB);
    #endif

    int nz = 0;
    for (int i = 0; i < n; i++) {
        // printf("Processing row %d\n", i);
        result->PI[i] = nz;
        for (int j = A->PI[i]; j < A->PI[i + 1]; j++) {
            result->J[nz] = A->J[j];
            result->VAL[nz] = scalarA * A->VAL[j];
            // printf("  Added A element: row %d, col %d, value = %f\n", i, A->J[j], result->VAL[nz]);
            nz++;
        }
        for (int j = B->PI[i]; j < B->PI[i + 1]; j++) {
            result->J[nz] = B->J[j];
            result->VAL[nz] = scalarB * B->VAL[j];
            // printf("  Added B element: row %d, col %d, value = %f\n", i, B->J[j], result->VAL[nz]);
            nz++;
        }
    }
    result->PI[n] = nz;
    result->nz = nz;

    #ifdef DEBUG_print_AddScaled
        // Final result structure
        printf("Final result matrix structure:\n");
        printf("PI array: ");
        for (int i = 0; i <= n; i++) {
            printf("%d ", result->PI[i]);
        }
        printf("\nJ and VAL arrays:\n");
        for (int i = 0; i < result->nz; i++) {
            printf("J[%d] = %d, VAL[%d] = %f\n", i, result->J[i], i, result->VAL[i]);
        }
        printf("CSR_AddScaled complete. Total non-zero entries: %d\n", result->nz);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_Consolidate(csr *A) {
    // Temporary storage for consolidated values
    int *temp_col = malloc(A->PI[A->n] * sizeof(int));
    double *temp_val = malloc(A->PI[A->n] * sizeof(double));
    int *temp_PI = malloc((A->n + 1) * sizeof(int));

    if (!temp_col || !temp_val || !temp_PI) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(temp_col);
        free(temp_val);
        free(temp_PI);
        return;
    }

    int nnz = 0;  // Count of non-zero entries after consolidation

    // Loop over each row in the matrix A
    for (int i = 0; i < A->n; i++) {
        temp_PI[i] = nnz;  // Start of row `i` in the consolidated matrix

        for (int j = A->PI[i]; j < A->PI[i + 1]; j++) {
            int col = A->J[j];
            double value = A->VAL[j];

            int found = 0;
            for (int k = temp_PI[i]; k < nnz; k++) {
                if (temp_col[k] == col) {
                    temp_val[k] += value;  // Consolidate values if the column index matches
                    found = 1;
                    break;
                }
            }

            // If column was not found in this row, add a new entry
            if (!found) {
                temp_col[nnz] = col;
                temp_val[nnz] = value;
                nnz++;
            }
        }
    }
    temp_PI[A->n] = nnz;  // End of the last row

    // Free original arrays and replace them with consolidated arrays
    free(A->J);
    free(A->VAL);
    free(A->PI);

    A->PI = malloc((A->n + 1) * sizeof(int));
    A->J = malloc(nnz * sizeof(int));
    A->VAL = malloc(nnz * sizeof(double));

    if (!A->PI || !A->J || !A->VAL) {
        fprintf(stderr, "Memory allocation failed for consolidated matrix.\n");
        free(temp_col);
        free(temp_val);
        free(temp_PI);
        return;
    }

    // Copy temporary arrays into the CSR structure of A
    memcpy(A->PI, temp_PI, (A->n + 1) * sizeof(int));
    memcpy(A->J, temp_col, nnz * sizeof(int));
    memcpy(A->VAL, temp_val, nnz * sizeof(double));

    // Free temporary memory
    free(temp_col);
    free(temp_val);
    free(temp_PI);
}
/* ------------------------------------------------------------------------------------------------- */
double CSR_RowResiduum(csr* A, const double *b, const double *x, int i, double *aii)
{
    #ifdef DEBUG_print_RowResiduum
    printf("Entering CSR_RowResiduum for row %d\n", i);
    printf("A->PI[%d] = %d, A->PI[%d + 1] = %d, A->nalloc = %d\n", i, A->PI[i], i, A->PI[i + 1], A->nalloc);
    #endif  
    
    // Verify that A->PI[i] and A->PI[i + 1] are within bounds
    if (i < 0 || i >= A->n || A->PI[i] < 0 || A->PI[i+1] > A->nalloc) {
        printf("Error: Array index out of bounds\n");
        // Consider adding error handling or exiting to avoid the segmentation fault
        return -1; // or appropriate error handling
    }

    double sum = 0.0;
    *aii = 0.0;

    // Loop through the non-zero elements of row 'i'
    for (int k = A->PI[i]; k < A->PI[i + 1]; k++)
    {
        #ifdef DEBUG_print_RowResiduum
        printf("Processing A->J[%d] within bounds (0, %d)\n", k, A->nalloc - 1);
        printf("k=%d, J[k]=%d, VAL[k]=%f, x[J[k]]=%f, Current sum=%f\n", k, A->J[k], A->VAL[k], x[A->J[k]], sum);
        if (A->J[k] < 0 || A->J[k] >= A->n) {
            printf("Error: A->J[k] is out of bounds\n");
            return -1; // or appropriate error handling
        }
        #endif

        if (A->J[k] == i)
        {
            *aii = A->VAL[k];  // Diagonal element
            sum += (*aii) * x[i];
            //printf("Sum[%d] = %g\n",k,sum);
        }
        else
        {
            sum += A->VAL[k] * x[A->J[k]];  // non-diagonal elements
        }
    }
    
    // printf("%d\n",i);

    #ifdef DEBUG_print_RowResiduum
    printf("Exiting CSR_RowResiduum for row %d, Residual: %f\n", i, b[i] - sum);
    #endif

    return b[i] - sum;  // Residual for row 'i'
}
/* ------------------------------------------------------------------------------------------------- */
void Gauss_Seidel(csr *A, const double *b, double *x, int max_iter, double tolerance)
{   
    #ifdef DEBUG_print_all
    printf("Entering Gauss_Seidel\n");
    #endif

    #ifdef TIMER
    clock_t start, end;
    double cpu_time_used = 0.0; 

    start = clock();
    #endif

    int iter;
    double final_rez = 0.0;
    int stopped_due_to_tolerance = 0;
    for (iter = 0; iter < max_iter; iter++)
    {
        double rezmax = 0;
        for (int i = 0; i < A->n; i++)
        {
            #ifdef DEBUG_print_all
            printf("Gauss_Seidel: Iteration %d\n", iter);
            #endif

            double aii = 0;
            double rezi = CSR_RowResiduum(A, b, x, i, &aii);

            if (aii != 0)
            {
                // Update x[i] using the residual and diagonal element
                double old_xi = x[i];
                x[i] = x[i] + rezi / aii;
                rezmax = fmax(rezmax, fabs(x[i] - old_xi));
            }
        }

        if (rezmax < tolerance)
        {
            stopped_due_to_tolerance = 1;
            final_rez = rezmax;
            #ifdef DEBUG_print_all
            printf("Tolerance reached at iteration %d\n", iter);
            #endif
            break;
        }
    }
    if (stopped_due_to_tolerance) {
        printf("Gauss-Seidel stopped due to reaching tolerance ( rez = %g ) ( iterations = %d ).\n", final_rez, iter);
    } else {
        printf("Gauss-Seidel stopped after reaching maximum iterations ( iter = %d ).\n",max_iter);
    }
    #ifdef DEBUG_print_all
    printf("Exiting Gauss_Seidel\n");
    #endif

    #ifdef TIMER
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time taken for Gauss-Seidel: %f seconds\n", cpu_time_used);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void SOR(csr *A, const double *b, double *x, int max_iter, double tolerance, double omega)
{
    #ifdef DEBUG_print_SOR
    printf("Entering SOR\n");
    #endif

    #ifdef TIMER
    clock_t start, end;
    double cpu_time_used = 0.0;

    start = clock();
    #endif

    //SOR - successive over-relaxation
    int iter;
    double final_rez = 0.0;
    int stopped_due_to_tolerance = 0;

    for (iter = 0; iter < max_iter; iter++)
    {
        #ifdef DEBUG_print_SOR
        printf("SOR: Iteration %d\n", iter);
        #endif
        double rezmax = 0;
        for (int i = 0; i < A->n; i++)
        {
            double aii = 0;
            double rezi = CSR_RowResiduum(A, b, x, i, &aii);

            // Check if diagonal element is non-zero
            if (aii != 0)
            {
                // Relaxation factor added
                double old_xi = x[i];
                x[i] = x[i] + (omega * rezi / aii);
                rezmax = fmax(rezmax, fabs(x[i] - old_xi));
                #ifdef DEBUG_print_SOR
                printf("x[%d]: %f\n",i,x[i]);
                printf("rezmax: %f\n",rezmax);
                printf("MAX_residuum: %f\n",rezmax);
                #endif
            }
        }
#ifdef PRNT_STEPS        
        if (iter % 500 == 0) {
            printf("Iteration %d, Residuum: %g\n", iter, rezmax);
        }
#endif
        if (rezmax < tolerance)
        {
            stopped_due_to_tolerance = 1;
            final_rez = rezmax;
            #ifdef DEBUG_print_SOR
            printf("Tolerance reached at iteration %d\n", iter);
            #endif
            break;
        }
    }
    if (stopped_due_to_tolerance) {
        printf("SOR stopped due to reaching tolerance ( rez = %g ) ( iterations = %d ).\n", final_rez, iter);
    } else {
        printf("SOR stopped after reaching maximum iterations ( iter = %d ).\n",max_iter);
    }

    #ifdef TIMER
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time taken for SOR: %f seconds\n", cpu_time_used);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void ConjugateGradient(csr* A, const double* b, double* x, int max_iter, double tolerance)
{
    #ifdef DEBUG_print_all
    printf("Entering ConjugateGradient\n");
    #endif

    #ifdef TIMER
    clock_t start, end;
    double cpu_time_used = 0.0;

    start = clock();
    #endif

    int n = A->n;
    double *r = malloc(n * sizeof(double));
    double *p = malloc(n * sizeof(double));
    double *Ap = malloc(n * sizeof(double));
    if (r == NULL || p == NULL || Ap == NULL) {
        printf("Memory Error ConjugateGradient: allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Initial residual r = b - Ax
    CSR_Multiply(A, x, Ap); // Computes Ap = A*x
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - Ap[i];
        p[i] = r[i];
    }

    double rsold = 0.0;
    for (int i = 0; i < n; i++) {
        rsold += r[i] * r[i];
    }

    int iter;
    int stopped_due_to_tolerance = 0;
    double final_rez = 0.0;

    for (iter = 0; iter < max_iter; iter++) {
        CSR_Multiply(A, p, Ap); // Ap = A*p
        double alpha = rsold;
        double Ap_dot_p = 0.0;
        for (int i = 0; i < n; i++) {
            Ap_dot_p += Ap[i] * p[i];
        }
        
        if (fabs(Ap_dot_p) < 1e-12) {
            fprintf(stderr, "Conjugate Gradient failed: Division by near-zero value Ap_dot_p.\n");
            break;
        }

        alpha /= Ap_dot_p;

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = 0.0;
        for (int i = 0; i < n; i++) {
            rsnew += r[i] * r[i];
        }
#ifdef PRNT_STEPS
        if (iter % 500 == 0) { // Add this condition to print every 20th iteration
            printf("Iteration %d, Residuum: %f\n", iter, sqrt(rsnew));
        }
#endif
        if (sqrt(rsnew) < tolerance) {
            stopped_due_to_tolerance = 1;
            final_rez = sqrt(rsnew);
            #ifdef DEBUG_print_all
            printf("ConjugateGradient reached tolerance at iteration %d\n", iter);
            #endif
            break;
        }

        for (int i = 0; i < n; i++) {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }

        rsold = rsnew;
    }

    if (!stopped_due_to_tolerance) {
        printf("ConjugateGradient stopped after reaching maximum iterations ( rez = %g ) ( iter = %d ).\n", final_rez, max_iter);
    } else {
        printf("ConjugateGradient stopped due to reaching tolerance ( rez = %g ) ( iterations = %d ).\n", final_rez, iter);
    }

    free(r);
    free(p);
    free(Ap);

    #ifdef TIMER
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time taken for Conjugate Gradient: %.10f seconds\n", cpu_time_used);
    #endif

    #ifdef DEBUG_print_all
    printf("Exiting ConjugateGradient\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void CSR_Solve(csr *A, const double *b, double *x) {
    if (!A || !b || !x) {
        fprintf(stderr, "Error: Null pointer in CSR_Solve inputs.\n");
        exit(EXIT_FAILURE);
    }

    // Convert CSR matrix for UMFPACK by treating it as transposed CSC
    void *Symbolic, *Numeric;

    // Perform symbolic factorization
    int status = umfpack_di_symbolic(A->n, A->n, A->PI, A->J, A->VAL, &Symbolic, NULL, NULL);
    if (status != UMFPACK_OK) {
        umfpack_di_report_status(NULL, status);
        fprintf(stderr, "UMFPACK symbolic factorization failed with status %d.\n", status);
        exit(EXIT_FAILURE);
    }

    // Perform numeric factorization
    status = umfpack_di_numeric(A->PI, A->J, A->VAL, Symbolic, &Numeric, NULL, NULL);
    if (status != UMFPACK_OK) {
        umfpack_di_report_status(NULL, status);
        umfpack_di_free_symbolic(&Symbolic);
        fprintf(stderr, "UMFPACK numeric factorization failed with status %d.\n", status);
        exit(EXIT_FAILURE);
    }
    umfpack_di_free_symbolic(&Symbolic);

    // Solve the system A * x = b using the numeric factorization
    status = umfpack_di_solve(UMFPACK_A, A->PI, A->J, A->VAL, x, b, Numeric, NULL, NULL);
    if (status != UMFPACK_OK) {
        umfpack_di_report_status(NULL, status);
        fprintf(stderr, "UMFPACK solve failed with status %d.\n", status);
        exit(EXIT_FAILURE);
    }
    umfpack_di_free_numeric(&Numeric);

    // printf("Solution completed using UMFPACK.\n");
}
/* ------------------------------------------------------------------------------------------------- */
void scaleMatricesForDirichlet(mesh *this, csr *M, csr *D, double omega, double eps1, double eps2, int nDOF_per_node, int useDamping) {
    double xi = (eps1 / (2 * omega)) + (eps2 * omega / 2);
    double nNodes = this->NPoints;

    printf("Scaling Dirichlet nodes with parameters:\n");
    printf("Omega: %f\n", omega);
    printf("Xi: %f\n", xi);
    printf("Mass scaling factor: %f\n", 1.0 / (omega * omega));
    printf("Damping scaling factor: %f\n", 2 * xi / omega);

    if (nDOF_per_node == 2) {
        for (int i = 0; i < nNodes; i++) {
            if (this->PointMark[i] == DIRICHLET_BC) {
                // Iterate over each DOF associated with the current node
                for (int dof = 0; dof < nDOF_per_node; dof++) {
                    int globalDOF = i * nDOF_per_node + dof;

                    // Adjust D matrix (damping matrix) if it exists
                    if (useDamping && D) {
                        for (int j = D->PI[globalDOF]; j < D->PI[globalDOF + 1]; j++) {
                            if (D->J[j] == globalDOF) {
                                D->VAL[j] = 2 * xi / omega;
                            }
                        }
                    }

                    // Adjust M matrix (mass matrix)
                    for (int j = M->PI[globalDOF]; j < M->PI[globalDOF + 1]; j++) {
                        if (M->J[j] == globalDOF) {
                            M->VAL[j] *= (1.0 / (omega * omega));
                        }
                    }
                }
            }
        }
    }
    else if (nDOF_per_node == 3) {
        for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
            if (this->ElementMark[elem] == DIRICHLET_BC) {
                gelement3D face;
                GetElement3D(this, &face, elem);
                int nVertices = GetNumBasisFunctions(face.type);
                
                for (int v = 0; v < nVertices; v++) {
                    int node = face.idxNode[v];
                    
                    for (int dof = 0; dof < nDOF_per_node; dof++) {
                        int globalDOF = node * nDOF_per_node + dof;

                        if (useDamping && D) {
                            for (int j = D->PI[globalDOF]; j < D->PI[globalDOF + 1]; j++) {
                                if (D->J[j] == globalDOF) {
                                    D->VAL[j] = 2 * xi / omega;
                                }
                            }
                        }
                        for (int j = M->PI[globalDOF]; j < M->PI[globalDOF + 1]; j++) {
                            if (M->J[j] == globalDOF) {
                                M->VAL[j] *= (1.0 / (omega * omega));
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        printf("Error: Unsupported nDOF_per_node = %d\n", nDOF_per_node);
        exit(EXIT_FAILURE);
    }
    
}

