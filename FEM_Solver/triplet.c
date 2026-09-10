#include "triplet.h"
#include "mesh.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define DEBUG_memory
// #define DEBUG_add
// #define DEBUG_sort
// #define DEBUG_unique
// #define DEBUG_finalize

void Triplet_Allocate(triplet* this, int nbNodes ,int nmax)
{
 
    #ifndef DEBUG_memory
    this->tripletI = malloc(nmax * sizeof(int));
    this->J = malloc(nmax * sizeof(int));
    this->VAL = malloc(nmax * sizeof(double));
    #endif

    #ifdef DEBUG_memory
    this->I = calloc(nmax , sizeof(int));
    this->J = calloc(nmax , sizeof(int));
    this->VAL = calloc(nmax , sizeof(double));
    #endif

    this->n = nbNodes;
    this->nz = 0;
    this->nalloc = nmax;

    #ifdef DEBUG_memory
    printf("Initial values after allocation:\n");
    for (int p = 0; p < nmax; p++)
    {
        printf("I: %d, J: %d, VAL: %f\n", this->I[p], this->J[p], this->VAL[p]);
    }
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Free(triplet* this)
{
    free(this->tripletI);
    free(this->J);
    free(this->VAL);

    this->nz = 0;
    this->n = 0;
    this->nalloc = 0;
}
/* ------------------------------------------------------------------------------------------------- */
void SaveTripletToFile(const char *filename, triplet *this) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Failed to open file %s for writing.\n", filename);
        return;
    }

    for (int i = 0; i < this->nz; i++) {
        fprintf(file, "%d %d %.15e\n", this->tripletI[i], this->J[i], this->VAL[i]);
    }

    fclose(file);
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Add(triplet* this, int i, int j, double aij)
{
 
    // printf("Add this->nz: %d\n", this->nz);

    if (this->nz >= this->nalloc) {
        fprintf(stderr, "Error: Attempting to add more elements than allocated space allows (Triplet_Add).\n");
        exit(EXIT_FAILURE);
        return;
    }

    #ifdef DEBUG_add
    printf("Soucet: I: %d, J: %d, VAL: %f\n", i, j, aij);
    #endif
    #ifdef PART
    if (this->nz >= this->nalloc) {
        int new_size = this->nalloc * 2;
        Triplet_Resize(this, new_size);
    }

    // Check for existing entries to avoid duplicates - ADDED THIS FOR LOOP
    for (int p = 0; p < this->nz; p++) {
        if (this->tripletI[p] == i && this->J[p] == j) {
            this->VAL[p] += aij;  // If the same i,j pair exists, add the values
            return;
        }
    }
    #endif

    this->tripletI[this->nz] = i;
    this->J[this->nz] = j;
    this->VAL[this->nz] = aij;

    this->nz++;

    #ifdef DEBUG_add
    printf("State after adding:\n");
    for (int p = 0; p < this->nz; p++)
    {
        printf("I: %d, J: %d, VAL: %f\n", this->I[p], this->J[p], this->VAL[p]);
    }
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Damping(triplet* D, const triplet* M, const triplet* K, const double eps1, const double eps2)
{
    
    Triplet_Allocate(D, M->n, M->nz + K->nz);

    for (int i = 0; i < M->nz; i++) {
        int row = M->tripletI[i];
        int col = M->J[i];
        double value = eps1 * M->VAL[i];
        Triplet_Add(D, row, col, value);
    }
    
    for (int i = 0; i < K->nz; i++) {
        int row = K->tripletI[i];
        int col = K->J[i];
        double value = eps2 * K->VAL[i];
        Triplet_Add(D, row, col, value);
    }

    Triplet_Finalize(D);
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Swap(triplet* this, int a, int b) 
{
    int tempI = this->tripletI[a];
    this->tripletI[a] = this->tripletI[b];
    this->tripletI[b] = tempI;

    int tempJ = this->J[a];
    this->J[a] = this->J[b];
    this->J[b] = tempJ;

    double tempVal = this->VAL[a];
    this->VAL[a] = this->VAL[b];
    this->VAL[b] = tempVal;
}
/* ------------------------------------------------------------------------------------------------- */
int Triplet_Partition(triplet* this, int first, int last)
{
    
    #ifdef DEBUG_sort
    printf("Triplet_Partition: first = %d, last = %d\n", first, last);
    #endif

    int pivotIndex = last;

    int pivotI = this->tripletI[pivotIndex];
    int pivotJ = this->J[pivotIndex];
    int i = first - 1;

    for (int j = first; j < last; j++)
    {
        if (this->tripletI[j] < pivotI || (this->tripletI[j] == pivotI && this->J[j] <= pivotJ)) //ADJUSTMENT
        {
            
            i++;

            #ifdef DEBUG_sort
            printf("Swapping: I[%d] = %d, J[%d] = %d with I[%d] = %d, J[%d] = %d\n", 
                   i, this->I[i], i, this->J[i], j, this->I[j], j, this->J[j]);
            #endif
            
            Triplet_Swap(this, i, j);
        }
    }

    #ifdef DEBUG_sort
    printf("Final swap before partitioning: I[%d] = %d, J[%d] = %d with I[%d] = %d, J[%d] = %d\n", 
           i + 1, this->I[i + 1], i + 1, this->J[i + 1], last, this->I[last], last, this->J[last]);
    #endif

    Triplet_Swap(this, i + 1, last);
    return i + 1;
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Sort(triplet* this, int first, int last)
{

    #ifdef DEBUG_sort
    printf("Triplet_Sort: Sorting range [%d, %d]\n", first, last);
    #endif

    if (first < last)
    {
    
        int pivotIndex = Triplet_Partition(this, first, last);

        #ifdef DEBUG_sort
        printf("Triplet_Sort: Pivot found at %d\n", pivotIndex);
        #endif

        if (pivotIndex > first)
        {
            #ifdef DEBUG_sort
            printf("Triplet_Sort: Sorting left half [%d, %d]\n", first, pivotIndex - 1);
            #endif

            Triplet_Sort(this, first, pivotIndex - 1);
        }
        
        if (pivotIndex < last)
        {
            #ifdef DEBUG_sort
            printf("Triplet_Sort: Sorting right half [%d, %d]\n", pivotIndex + 1, last);
            #endif

            Triplet_Sort(this, pivotIndex + 1, last);
        }
    }

    // printf("Triplet_Sort: tri.nz: %d \t last: %d\n", this->nz, last);

    #ifdef DEBUG_sort
    printf("Triplet_Sort: Sorted range [%d, %d]\n", first, last);
    for (int p = first; p <= last; p++)
    {
        printf("I: %d, J: %d, VAL: %f\n", this->tripletI[p], this->J[p], this->VAL[p]);
    }
    #endif

    // printf("After sorting:\n");
    // for (int p = first; p <= last; p++) {
    //     printf("I: %d, J: %d, VAL: %f\n", this->tripletI[p], this->J[p], this->VAL[p]);
    // }
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Sort1(triplet* this, int first, int last)
{
    int nz = this->nz;
    triplet_entry* temp = malloc(nz * sizeof(triplet_entry));
    if (!temp) {
        fprintf(stderr, "Triplet_Sort: malloc failed.\n");
        exit(EXIT_FAILURE);
    }

    for (int k = 0; k < nz; ++k) {
        temp[k].i = this->tripletI[k];
        temp[k].j = this->J[k];
        temp[k].val = this->VAL[k];
    }

    qsort(temp, nz, sizeof(triplet_entry), compare_triplet_entries);

    for (int k = 0; k < nz; ++k) {
        this->tripletI[k] = temp[k].i;
        this->J[k] = temp[k].j;
        this->VAL[k] = temp[k].val;
    }

    free(temp);
}
/* ------------------------------------------------------------------------------------------------- */
int compare_triplet_entries(const void* a, const void* b)
{
    const triplet_entry* ta = (const triplet_entry*)a;
    const triplet_entry* tb = (const triplet_entry*)b;

    if (ta->i < tb->i) return -1;
    if (ta->i > tb->i) return 1;
    if (ta->j < tb->j) return -1;
    if (ta->j > tb->j) return 1;
    return 0;
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Unique(triplet* this) 
{
    #ifdef DEBUG_unique
    printf("Before merging:\n");
    for (int p = 0; p < this->nz; p++)
    {
        printf("I: %d, J: %d, VAL: %f\n", this->tripletI[p], this->J[p], this->VAL[p]);
    }
    #endif

    // SaveTripletToFile("Triplet_Before_Unique_full.dat", this);

    int i = 0, j = 1;
    for ( i = 0; i < this->nz; i++) {
    //    while (this->tripletI[i] == this->tripletI[j] && this->J[i] == this->J[j]) {
    //         this->VAL[i] += this->VAL[j];
    //         j++;
    //     }
        while (j < this->nz && this->tripletI[i] == this->tripletI[j] && this->J[i] == this->J[j]) {
            this->VAL[i] += this->VAL[j];
            j++;
        }        
        if (j < this->nz) {
            this->tripletI[i + 1]   = this->tripletI[j];
            this->J[i + 1]   = this->J[j];
            this->VAL[i + 1] = this->VAL[j];
            j++;
        }
        else
            break;
    }
    this->nz = i + 1;
    
    // SaveTripletToFile("Triplet_After_Unique_full.dat", this);

    // #define DEBUG_unique
    #ifdef DEBUG_unique
    printf("After merging of Triplet_Unique:\n");
    for (int p = 0; p < this->nz; p++)
    {
        printf("I: %d, J: %d, VAL: %f\n", this->tripletI[p], this->J[p], this->VAL[p]);
    }
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Set(triplet* this, int i, int j, double aij)
{
    
    int found = 0;

    // Search for existing entry with the same indices
    for (int p = 0; p < this->nz; p++) {
        if (this->tripletI[p] == i && this->J[p] == j) {
            this->VAL[p] = aij; // Update the existing entry
            found = 1;
            break;
        }
    }

    // If no existing entry was found, add a new one
    if (!found) {
        Triplet_Add(this, i, j, aij);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Resize(triplet* this, int new_size) {
    //*
    if (new_size <= this->nalloc) {
        return;  // No need to resize if the new size is not larger
    }
    //*/
    this->tripletI = realloc(this->tripletI, new_size * sizeof(int));
    this->J = realloc(this->J, new_size * sizeof(int));
    this->VAL = realloc(this->VAL, new_size * sizeof(double));
    
    if (!this->tripletI || !this->J || !this->VAL) {
        fprintf(stderr, "Error reallocating triplet arrays.\n");
        exit(EXIT_FAILURE);
    }
    
    this->nalloc = new_size;
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Finalize(triplet* this)
{
    if (this->nz > 0)
    {
        Triplet_Sort(this, 0, this->nz - 1);
        Triplet_Unique(this);

    }
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Reset(triplet* this)
{
    this->nz = 0;
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Save(const triplet* this, const char* fname) {
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Triplet_Save]: Failed to open file for saving TRIPLET data.\n");
        exit(1);
    }

    for (int i = 0; i < this->nz; i++) {
        fprintf(fid, "%d %d %e\n", this->tripletI[i] + 1, this->J[i] + 1, this->VAL[i]);
    }

    fclose(fid);
    printf("Data saved successfully to %s with 1-based indexing.\n", fname);
}
/* ------------------------------------------------------------------------------------------------- */
void Triplet_Init(triplet* this, int n, int m, int nmax) {
    Triplet_Allocate(this, n, nmax);  
    this->nz = 0;                     
    this->n = n;                      
}
/* ------------------------------------------------------------------------------------------------- */


