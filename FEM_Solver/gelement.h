#ifndef GELEMENT_H
#define GELEMENT_H

#include "mesh.h"

typedef struct{

    int iA;
    int iB;
    int iC;
    int iD;
    
    double A[2];
    double B[2];
    double C[2];
    double D[2];

    double Bk[2][2];

    double detBk;

    double invBk[2][2];

    int idx[3];
    int idxEle;
    double dV, dS;
    
}gelement;

gelement Gelement(const mesh *M, int i);

double calculateVolume(const mesh *M);

//double calculateTriangleArea(double A[2], double B[2], double C[2]);

#endif
