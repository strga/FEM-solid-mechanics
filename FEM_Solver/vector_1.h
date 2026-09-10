#pragma once

typedef double* vector;

vector Vector_Allocate(int n);
void   Vector_Free(vector x);
void   Vector_Linspace(vector* x, int start, double dx, int end, int* n);
void   Vector_Save(vector x, int n, const char* fname);
void   Vector_Save2(vector x, vector y, int n, const char* fname);
void   Vector_Save3(vector x, vector y, vector z, int n, const char* fname);
void   Vector_Save4(vector x, vector y, vector z, vector z_, int n, const char* fname);

void   Vector_EuclidNorm(vector x, double* Enorm, int n);
void   Vector_EuclidINorm(vector a, vector b, double* Enorm, int n);

void   Matrix_Save(double arr[][6], int n, const char *fname);
void   Matrix_Save2(double arr[][12], int n, const char *fname);
void   Matrix_Save3(double arr[][7], int n, const char *fname);
