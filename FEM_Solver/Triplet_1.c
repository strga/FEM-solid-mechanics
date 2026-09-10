//#define DEBUG
#include <stdlib.h>
#include <stdio.h>
#include "Triplet_1.h"


void Triplet_1_Allocate(triplet_1* this, int nbNodes, int nmax) {
    this->I_re = malloc(nmax * sizeof(int));
    if (this->I_re == NULL) {
        printf("\nError[Triplet_Allocate]: Nedostatek pameti pro triplet.I");
        exit(1);
    }

    this->J = malloc(nmax * sizeof(int));
    if (this->J == NULL) {
        printf("\nError[Triplet_Allocate]: Nedostatek pameti pro triplet.J");
        exit(1);
    }

    this->VAL = malloc(nmax * sizeof(double));
    if (this->VAL == NULL) {
        printf("\nError[Triplet_Allocate]: Nedostatek pameti pro triplet.VAL");
        exit(1);
    }

    this->n      = nbNodes;
    this->nz     = 0;
    this->nalloc = nmax;
}
/* ----------------------------------------------------------------------------------- */
void Triplet_1_Free(triplet_1* this) {
    if (this->nalloc > 0) {
        free(this->I_re);
        free(this->J);
        free(this->VAL);

        this->nalloc = 0;
        this->nz     = 0;
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet_1_Add(triplet_1* this, int i, int j, double val) {
    
    // printf("Add this1->nz: %d\n", this->nz);
    
    if (this->nz >= this->nalloc) {
        printf("\nError[Triplet_Add]: Alokovano malo pameti pro TRIPLET\n");
        exit(1);
    }
    this->I_re[this->nz]   = i;
    this->J[this->nz]   = j;
    this->VAL[this->nz] = val;

#ifdef DEBUG
    printf("I = %d J = %d VAL = %lf\n", this->I[this->nz], this->J[this->nz], this->VAL[this->nz]);
#endif // DEBUG

    this->nz++;
}
/* ----------------------------------------------------------------------------------- */
/*
void Triplet_1_AddRHS(triplet_1* this, vector b, int i, int j, double qval) {
    if (i > 0) {
        if (j > 0) {
            Triplet_Add(this, i, j, qval);
        }
        else if (j < 0) {
            b[i - 1] -= qval * b[-j - 1];
        }
    }
}
*/
/* ----------------------------------------------------------------------------------- */
void Triplet_1_Quicksort(triplet_1* this, int first, int last) {
    int i, pivot;
    double temp;

    if (first < last) {
        // printf("first = %d/ last = %d\n", first, last);
        pivot = first;
        for (i = first + 1; i < last; i++) {
            if (this->I_re[i] < this->I_re[first]) {
                pivot++;

                temp = this->I_re[pivot];
                this->I_re[pivot] = this->I_re[i];
                this->I_re[i] = temp;

                temp = this->J[pivot];
                this->J[pivot] = this->J[i];
                this->J[i] = temp;

                temp = this->VAL[pivot];
                this->VAL[pivot] = this->VAL[i];
                this->VAL[i] = temp;
            }
            if (this->I_re[i] == this->I_re[first] && this->J[i] < this->J[first]) {
                pivot++;

                temp = this->I_re[pivot];
                this->I_re[pivot] = this->I_re[i];
                this->I_re[i] = temp;

                temp = this->J[pivot];
                this->J[pivot] = this->J[i];
                this->J[i] = temp;

                temp = this->VAL[pivot];
                this->VAL[pivot] = this->VAL[i];
                this->VAL[i] = temp;
            }
        }

        temp = this->I_re[pivot];
        this->I_re[pivot] = this->I_re[first];
        this->I_re[first] = temp;

        temp = this->J[pivot];
        this->J[pivot] = this->J[first];
        this->J[first] = temp;

        temp = this->VAL[pivot];
        this->VAL[pivot] = this->VAL[first];
        this->VAL[first] = temp;

        Triplet_1_Quicksort(this, first, pivot);
        Triplet_1_Quicksort(this, pivot + 1, last);
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet_1_Unique(triplet_1* this) {
    int i, j = 1;
    printf("\n");
    printf("Triplet_1_Unique BEFORE: this.nz: %d\n", this->nz);

    // printf("Values of I_re[]:\n");
    // for (int i = 0; i < this->nz; i++) {
    //     printf("I_re[%d] = %d\n", i, this->I_re[i]);
    // }

    for (i = 0; i < this->nz; i++) {
        while (this->I_re[i] == this->I_re[j] && this->J[i] == this->J[j]) {
            this->VAL[i] += this->VAL[j];
            j++;
        }
        if (j < this->nz) {
            this->I_re[i + 1]   = this->I_re[j];
            this->J[i + 1]   = this->J[j];
            this->VAL[i + 1] = this->VAL[j];
            j++;
        }
        else
            break;
    }
    this->nz = i + 1;
    printf("Triplet_1_Unique AFTER: this.nz:: %d\n",this->nz);
    
    // printf("\n");
    // printf("Values2 of I_re[]:\n");
    // for (int i = 0; i < this->nz; i++) {
    //     printf("I_re[%d] = %d\n", i, this->I_re[i]);
    // }

    // #define DEBUG_unique
    #ifdef DEBUG_unique
    printf("After merging of Triplet_Unique:\n");
    for (int p = 0; p < this->nz; p++)
    {
        printf("I: %d, J: %d, VAL: %f\n", this->I_re[p], this->J[p], this->VAL[p]);
    }
    #endif
}
/* ----------------------------------------------------------------------------------- */
void Triplet_1_Save(const triplet_1* this, const char* fname) {
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Triplet_Save]: Nepodarilo se otevrit soubor pro TRIPLET");
        exit(1);
    }

    for (int i = 0; i < this->nz; i++) 
        fprintf(fid, "%d %d %lf\n", this->I_re[i], this->J[i], this->VAL[i]);
        
    fclose(fid);
}


void Triplet_Imag_Allocate(itriplet* this, int nbNodes, int nmax) {
    this->I_im = malloc(nmax * sizeof(int));
    if (this->I_im == NULL) {
        printf("\nError[Triplet_Imag_Allocate]: Nedostatek pameti pro triplet.I");
        exit(1);
    }

    this->J = malloc(nmax * sizeof(int));
    if (this->J == NULL) {
        printf("\nError[Triplet_Imag_Allocate]: Nedostatek pameti pro triplet.J");
        exit(1);
    }

    this->VAL = malloc(nmax * sizeof(double));
    if (this->VAL == NULL) {
        printf("\nError[Triplet_Imag_Allocate]: Nedostatek pameti pro triplet.VAL");
        exit(1);
    }

    this->IVAL = malloc(nmax * sizeof(double));
    if (this->IVAL == NULL) {
        printf("\nError[Triplet_Imag_Allocate]: Nedostatek pameti pro triplet.IVAL");
        exit(1);
    }

    this->n      = nbNodes;
    this->nz     = 0;
    this->nalloc = nmax;
}
/* ----------------------------------------------------------------------------------- */
void Triplet_Imag_Free(itriplet* this) {
    if (this->nalloc > 0) {
        free(this->I_im);
        free(this->J);
        free(this->VAL);
        free(this->IVAL);

        this->nalloc = 0;
        this->nz     = 0;
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet_Imag_Add(itriplet* this, int i, int j, double val, double ival) {
    if (this->nz >= this->nalloc) {
        printf("\nError[Triplet_Imag_Add]: Alokovano malo pameti pro ITRIPLET\n");
        exit(1);
    }
    this->I_im[this->nz]    = i;
    this->J[this->nz]    = j;
    this->VAL[this->nz]  = val;
    this->IVAL[this->nz] = ival;

#ifdef DEBUG
    printf("I = %d J = %d VAL = %lf IVAL = %lf\n", this->I_im[this->nz], this->J[this->nz], this->VAL[this->nz], this->IVAL[this->nz]);
#endif // DEBUG

    this->nz++;
}
/* ----------------------------------------------------------------------------------- */
void Triplet_Imag_AddRHS(itriplet* this, vector b, vector ib, int i, int j, double qval, double iqval) {
    if (i > 0) {
        if (j > 0) {
            Triplet_Imag_Add(this, i, j, qval, iqval);
        }
        else if (j < 0) {
            b[i - 1]  -= qval * b[-j - 1];
            ib[i - 1] -= qval * ib[-j - 1];
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet_Imag_Quicksort(itriplet* this, int first, int last) {
    int i, pivot;
    double temp;

    if (first < last) {
        pivot = first;
        for (i = first + 1; i < last; i++) {
            if (this->I_im[i] < this->I_im[first]) {
                pivot++;

                temp = this->I_im[pivot];
                this->I_im[pivot] = this->I_im[i];
                this->I_im[i] = temp;

                temp = this->J[pivot];
                this->J[pivot] = this->J[i];
                this->J[i] = temp;

                temp = this->VAL[pivot];
                this->VAL[pivot] = this->VAL[i];
                this->VAL[i] = temp;

                temp = this->IVAL[pivot];
                this->IVAL[pivot] = this->IVAL[i];
                this->IVAL[i] = temp;
            }
            if (this->I_im[i] == this->I_im[first] && this->J[i] < this->J[first]) {
                pivot++;

                temp = this->I_im[pivot];
                this->I_im[pivot] = this->I_im[i];
                this->I_im[i] = temp;

                temp = this->J[pivot];
                this->J[pivot] = this->J[i]; 
                this->J[i] = temp;

                temp = this->VAL[pivot];
                this->VAL[pivot] = this->VAL[i];
                this->VAL[i] = temp;

                temp = this->IVAL[pivot];
                this->IVAL[pivot] = this->IVAL[i];
                this->IVAL[i] = temp;
            }
        }

        temp = this->I_im[pivot];
        this->I_im[pivot] = this->I_im[first];
        this->I_im[first] = temp;

        temp = this->J[pivot];
        this->J[pivot] = this->J[first];
        this->J[first] = temp;

        temp = this->VAL[pivot];
        this->VAL[pivot] = this->VAL[first];
        this->VAL[first] = temp;

        temp = this->IVAL[pivot];
        this->IVAL[pivot] = this->IVAL[first];
        this->IVAL[first] = temp;

        Triplet_Imag_Quicksort(this, first, pivot);
        Triplet_Imag_Quicksort(this, pivot + 1, last);
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet_Imag_Unique(itriplet* this) {
    int i, j = 1;

    for (i = 0; i < this->nz; i++) {
        while (this->I_im[i] == this->I_im[j] && this->J[i] == this->J[j]) {
            this->VAL[i]  += this->VAL[j];
            this->IVAL[i] += this->IVAL[j];
            j++;
        }
        if (j < this->nz) {
            this->I_im[i + 1]    = this->I_im[j];
            this->J[i + 1]    = this->J[j];
            this->VAL[i + 1]  = this->VAL[j];
            this->IVAL[i + 1] = this->IVAL[j];
            j++;
        }
        else
            break;
    }
    this->nz = i + 1;
}
/* ----------------------------------------------------------------------------------- */
void Triplet_Imag_Save(const itriplet* this, const char* fname) {
    FILE* fid;

    fid = fopen(fname, "w");
    if (fid == NULL) {
        printf("\nError[Triplet_Imag_Save]: Nepodarilo se otevrit soubor pro TRIPLET");
        exit(1);
    }

    for (int i = 0; i < this->nz; i++) 
        fprintf(fid, "%d %d %lf %lf\n", this->I_im[i], this->J[i], this->VAL[i], this->IVAL[i]);
        
    fclose(fid);
}
