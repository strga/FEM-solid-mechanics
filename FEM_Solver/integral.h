
#pragma once

//#ifndef INTEGRAL_H
//#define INTEGRAL_H

//#ifdef INTEGRATE_OMEGA_H
//#define INTEGRATE_OMEGA_H

#include "mesh.h"

//typedef double (*Integrate_Function)(double x, double y);

typedef struct
{
    double vertexValue;
    double midpointValue;
    double centroidValue;
    
} ElementIntegrationResult;


typedef struct
{
    
    double sumVertex;
    double sumMidPts;
    double sumTeziste;

} Integral;

typedef struct
{
    
    double sum1;
    double sum2;
    double sum3;

} Integrate_Omega;


Integral Mesh_Integrate(mesh *this, double f(double x, double y), int elementIdx);
Integrate_Omega Region_Integrate(mesh *this, double f(double x, double y));

//#endif