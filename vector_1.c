#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vector_1.h"


vector Vector_Allocate(int n) {
    vector p;

    p = malloc(n * sizeof(double));

    return p;
}
/* ----------------------------------------------------------------------------------- */
void Vector_Free(vector x) {
    free(x);
}
/* ----------------------------------------------------------------------------------- */
void Vector_Linspace(vector* x, int start, double dx, int end, int* n) {
    double n_ = (end - start) / dx;
    
    if (round(n_) != n_) {
        printf("\nError[Vector_Linspace]: Interval nelze rozdelit rovnomerne, pocet dilku: %lf\n", n_);
        exit(1);
    }

    *n = n_ + 1;
    *x = Vector_Allocate(*n);

    (*x)[0]    = start;
    (*x)[*n-1] = end;

    for(int i = 1; i < *n - 1; i++)
        (*x)[i] = (*x)[i-1] + dx; 
}
/* ----------------------------------------------------------------------------------- */
void Vector_Save(vector x, int n, const char* fname) {
    int   i;
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Vector_Save]: Nepodarilo se otevrit soubor pro VECTOR\n");
        exit(1);
    }

    for (i = 0; i < n; i++) {
        fprintf(fid, "%e\n", x[i]);
    }

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Vector_Save2(vector x, vector y, int n, const char* fname) {
    int   i;
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Vector_Save]: Nepodarilo se otevrit soubor pro VECTOR\n");
        exit(1);
    }

    for (i = 0; i < n; i++) {
        fprintf(fid, "%e %e\n", x[i], y[i]);
    }

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Vector_Save3(vector x, vector y, vector z, int n, const char* fname) {
    int   i;
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Vector_Save]: Nepodarilo se otevrit soubor pro VECTOR\n");
        exit(1);
    }

    for (i = 0; i < n; i++)
        fprintf(fid, "%e %e %e\n", x[i], y[i], z[i]);

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Vector_Save4(vector x, vector y, vector z, vector z_, int n, const char* fname) {
    int   i;
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Vector_Save]: Nepodarilo se otevrit soubor pro VECTOR\n");
        exit(1);
    }

    for (i = 0; i < n; i++)
        fprintf(fid, "%e %e %e %e\n", x[i], y[i], z[i], z_[i]);

    fclose(fid);
}



void Matrix_Save(double arr[][6], int n, const char *fname) {
    int i, j;
    FILE *fid;

    fid = fopen(fname, "w");
    if (fid == NULL){
        printf("\nError[Matrix_Save]: Nepodarilo se otevrit soubor\n");
        exit(1);
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < 6; j++)
            fprintf(fid, "%16.8e\t", arr[i][j]);
        fprintf(fid,"\n");
    }

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Matrix_Save2(double arr[][12], int n, const char *fname) {
    int i, j;
    FILE *fid;

    fid = fopen(fname, "w");
    if (fid == NULL){
        printf("\nError[Matrix_Save2]: Nepodarilo se otevrit soubor\n");
        exit(1);
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < 12; j++)
            fprintf(fid, "%16.8e\t", arr[i][j]);
        fprintf(fid,"\n");
    }

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Matrix_Save3(double arr[][7], int n, const char *fname) {
    int i, j;
    FILE *fid;

    fid = fopen(fname, "w");
    if (fid == NULL){
        printf("\nError[Matrix_Save3]: Nepodarilo se otevrit soubor\n");
        exit(1);
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < 7; j++)
            fprintf(fid, "%16.8e\t", arr[i][j]);
        fprintf(fid,"\n");
    }

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Vector_EuclidNorm(vector x, double* Enorm, int n) {
    *Enorm = 0;

    for (int i = 0; i < n; i++)
        *Enorm += pow(x[i], 2); 

    *Enorm = sqrt(*Enorm);
    //printf("Enorm = %e\n", *Enorm);
}
/* ----------------------------------------------------------------------------------- */
void Vector_EuclidINorm(vector a, vector b, double* Enorm, int n) {
    *Enorm = 0;

    for (int i = 0; i < n; i++)
        *Enorm += (pow(a[i], 2) + pow(b[i], 2)); 

    *Enorm = sqrt(*Enorm);
} 