#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "gelement.h"
#include "integral.h"

//#define PRINT

Integral Mesh_Integrate(mesh *this, double f(double x, double y), int elementIdx )
{
    
    Integral result;
    result.sumVertex = 0.0;
    result.sumMidPts = 0.0;
    result.sumTeziste = 0.0;

    gelement tri = Gelement(this, elementIdx);

    //double V1 = (f(tri.A[0], tri.A[1]) / 3.0 + f(tri.B[0], tri.B[1]) / 3.0  + f(tri.C[0], tri.C[1]) / 3.0 );
    double V1 = (f(tri.A[0], tri.A[1]) + f(tri.B[0], tri.B[1])  + f(tri.C[0], tri.C[1]) ) / 3.0;

    double MidA[] = {(tri.B[0] + tri.C[0]) / 2.0, (tri.B[1] + tri.C[1]) / 2.0};
    double MidB[] = {(tri.A[0] + tri.C[0]) / 2.0, (tri.A[1] + tri.C[1]) / 2.0};
    double MidC[] = {(tri.A[0] + tri.B[0]) / 2.0, (tri.A[1] + tri.B[1]) / 2.0};

    double V2 = (f(MidA[0], MidA[1])  + f(MidB[0], MidB[1])  + f(MidC[0], MidC[1]) ) / 3.0;

    double T[] = {(tri.A[0] + tri.B[0] + tri.C[0]) / 3.0, (tri.A[1] + tri.B[1] + tri.C[1])};

    double V3 = f(T[0], T[1]);

    result.sumVertex = V1 * tri.dV; // Vertex
    result.sumMidPts = V2 * tri.dV; // Mid points
    result.sumTeziste = V3 * tri.dV; // Centre of gravity

#ifdef PRINT
    printf("Element %d Integration Results:\n", elementIdx);
    printf("Area of element: %f\n",tri.dV);
    printf("  V1 (Vertices): %f\n", V1);
    printf("  V2 (Mid Points): %f\n", V2);
    printf("  V3 (Centroid): %f\n", V3);
    printf("\n");
#endif

    return result;
}

Integrate_Omega Region_Integrate(mesh *this, double f(double x, double y) )
{

    Integrate_Omega result;
    result.sum1 = 0.0;
    result.sum2 = 0.0;
    result.sum3 = 0.0;

    
    for (int i = 0; i < this->NElements; ++i) {
        gelement tri = Gelement(this, i);

        double W1 = (f(tri.A[0], tri.A[1]) + f(tri.B[0], tri.B[1]) + f(tri.C[0], tri.C[1])) / 3.0;

        double MidA[] = {(tri.B[0] + tri.C[0]) / 2.0, (tri.B[1] + tri.C[1]) / 2.0};
        double MidB[] = {(tri.A[0] + tri.C[0]) / 2.0, (tri.A[1] + tri.C[1]) / 2.0};
        double MidC[] = {(tri.A[0] + tri.B[0]) / 2.0, (tri.A[1] + tri.B[1]) / 2.0};
        
        double W2 = (f(MidA[0], MidA[1]) + f(MidB[0], MidB[1]) + f(MidC[0], MidC[1])) / 3.0;

        double T[] = {(tri.A[0] + tri.B[0] + tri.C[0]) / 3.0, (tri.A[1] + tri.B[1] + tri.C[1]) / 3.0};

        double W3 = f(T[0], T[1]);
        
        result.sum1 += W1 * tri.dV; // Vertex
        result.sum2 += W2 * tri.dV; // Mid points
        result.sum3 += W3 * tri.dV; // Centre of gravity
    }
    
    return result;
}

