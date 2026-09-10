#include "Mesh_GetBndrSide.h"
#include "mesh.h"
#include <math.h>

Mesh_GetBndrSide MeshGetBoundarySide(const mesh *M, int i) {
    Mesh_GetBndrSide S;

    S.idxA = M->EdgeVerA[i];
    S.idxB = M->EdgeVerB[i];
    S.mark = M->EdgeMark[i];

    double dx = M->x[S.idxB] - M->x[S.idxA];
    double dy = M->y[S.idxB] - M->y[S.idxA];

    S.ds = sqrt( (dx * dx) + (dy * dy) );

    S.nn[0] = dy / S.ds;
    S.nn[1] =  -dx / S.ds;

    return S;
}

double calculateBndrLength(const mesh *M)
{
    double Length = 0;

    for (int i = 0; i < M->NBEdges; i++)
    {
        Mesh_GetBndrSide S = MeshGetBoundarySide(M, i);
        Length += S.ds;
    }

    return Length;
}