#ifndef TRIPLET_H
#define TRIPLET_H

#include "mesh.h"

typedef struct{
    
    int *tripletI;
    int *J;
    double *VAL; // aij
    
    int nz;
    int n;
    int nalloc;

}triplet;

typedef struct {
    int i, j;
    double val;
} triplet_entry;


void Triplet_Allocate(triplet* this, int nbNodes,int nmax);
void Triplet_Add(triplet* this, int i, int j, double aij);

void Triplet_Sort1(triplet* this, int first, int last);
int compare_triplet_entries(const void* a, const void* b);

void Triplet_Sort(triplet* this, int first, int last);
void Triplet_Unique(triplet* this);
void Triplet_Set(triplet* this, int i, int j, double aij);
void Triplet_Finalize(triplet* this);
void Triplet_Save(const triplet* this, const char* fname);
void Triplet_Resize(triplet* this, int new_size);
void Triplet_Damping(triplet* D, const triplet* M, const triplet* K, const double eps1, const double eps2);
void Triplet_Free(triplet *this);
void Triplet_Swap(triplet* this, int a, int b);
int Triplet_Partition(triplet* this, int first, int last);
void Triplet_Reset(triplet* this);

void Triplet_Save(const triplet* this, const char* fname);
void PrintTriplet(const char *label, triplet *this);

void Triplet_Init(triplet* this, int n, int m, int nmax);
#endif