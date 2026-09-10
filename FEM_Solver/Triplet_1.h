#pragma once

#include "vector_1.h"

typedef struct {
    int      *I_re;        // indexy radku
    int      *J;        // indexy sloupcu
    double   *VAL;      // prvky matice
    int      n;         // rozmer matice
    int      nz;        // pocet nenulovych prvku
    int      nalloc;    // alokovana pamet
} triplet_1;

typedef struct {
    int      *I_im;        
    int      *J;        
    double   *VAL;      // realna cast
    double   *IVAL;     // imag cast
    int      n;         
    int      nz;        
    int      nalloc;    
} itriplet;

void Triplet_1_Allocate(triplet_1* this, int nbNodes, int nmax);
void Triplet_1_Free(triplet_1* this);
void Triplet_1_Add(triplet_1* this, int i, int j, double val);
// void Triplet_1_AddRHS(triplet_1* this, vector b, int i, int j, double qval);
void Triplet_1_Quicksort(triplet_1* this, int first, int last);
void Triplet_1_Unique(triplet_1* this);
void Triplet_1_Save(const triplet_1* this, const char* fname);


void Triplet_Imag_Allocate(itriplet* this, int nbNodes, int nmax);
void Triplet_Imag_Free(itriplet* this);
void Triplet_Imag_Add(itriplet* this, int i, int j, double val, double ival);
void Triplet_Imag_AddRHS(itriplet* this, vector b, vector ib, int i, int j, double qval, double iqval);
void Triplet_Imag_Quicksort(itriplet* this, int first, int last);
void Triplet_Imag_Unique(itriplet* this);
void Triplet_Imag_Save(const itriplet* this, const char* fname);