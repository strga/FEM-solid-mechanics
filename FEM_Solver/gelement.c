#include "gelement.h"
#include "mesh.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

gelement Gelement(const mesh *M, int i)
{
    
     if (i < 0 || i >= M->NElements) {
        fprintf(stderr, "Error: Element index %d is out of bounds.\n", i);
        exit(EXIT_FAILURE);
    }

    if (M->ElementListOfVertices == NULL || M->x == NULL || M->y == NULL) {
        fprintf(stderr, "Error: Null pointer encountered in mesh structure (gelement.c).\n");
        exit(EXIT_FAILURE);
    }

    gelement K;

    K.iA = M->ElementListOfVertices[3 * i + 0];
    K.iB = M->ElementListOfVertices[3 * i + 1];
    K.iC = M->ElementListOfVertices[3 * i + 2]; 

    K.idx[0] = K.iA;
    K.idx[1] = K.iB;
    K.idx[2] = K.iC;

    // printf("i:%d A:%d B:%d C:%d\n",i,K.iA,K.iB,K.iC);

    K.A[0] = M->x[K.iA];
    K.A[1] = M->y[K.iA];

    K.B[0] = M->x[K.iB];
    K.B[1] = M->y[K.iB];

    K.C[0] = M->x[K.iC];
    K.C[1] = M->y[K.iC];

    K.Bk[0][0] = K.B[0] - K.A[0];
    K.Bk[1][0] = K.B[1] - K.A[1];

    K.Bk[0][1] = K.C[0] - K.A[0];
    K.Bk[1][1] = K.C[1] - K.A[1];

    K.detBk = K.Bk[0][0] * K.Bk[1][1] - K.Bk[0][1] * K.Bk[1][0];

    K.dV = 0.5 * fabs(K.detBk);

    double det_inv = 1 / K.detBk;

    K.invBk[0][0] = K.Bk[1][1] * det_inv;
    K.invBk[0][1] = - K.Bk[0][1] * det_inv;

    K.invBk[1][0] = - K.Bk[1][0] * det_inv;
    K.invBk[1][1] = K.Bk[0][0] * det_inv;

    return K;
}

double calculateVolume(const mesh *M)
{
    double volume = 0;

    for (int i = 0; i < M->NElements; i++)
    {
        gelement K = Gelement(M, i);
        volume += K.dV;
    }
    
    return volume;
}