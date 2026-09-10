#ifndef MESHGETBNDRSIDE_H
#define MESHGETBNDRSIDE_H

#include "mesh.h"

typedef struct {
    
    int idxA, idxB, mark;
    double *A[2];
    double *B[2];
    double MidP[2];
    double dx;
    double dy;
    double ds;
    double nn[2];

} Mesh_GetBndrSide;

Mesh_GetBndrSide MeshGetBoundarySide(const mesh *M, int i);

double calculateBndrLength(const mesh *M);

#endif
