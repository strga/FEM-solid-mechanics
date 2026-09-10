#ifndef GELEMENT3D_H
#define GELEMENT3D_H

#include "mesh.h"

#define TRIANGLE_VERTEX      3  // TRIANGLE
#define QUADRANGLE_VERTEX    4  // QUADRANGLE
#define TETRAHEDRON_VERTEX   4  // TETRAHEDRON
#define HEXAHEDRON_VERTEX    8  // HEXAHEDRON
#define PYRAMID_VERTEX       5  // PYRAMID
#define PRISM_VERTEX         6  // PRISM


typedef struct {
    double x, y, z;
    int mark;
} vertex;

typedef struct {
    int type;               // Type of the geometric element
    int Nvertex;            // Number of vertices for given element
    int mark;               // Boundary condition marker
    int idxNode[8];         // Vertex indices
    vertex Ver[8];
    double volume;          // Volume of the element
    double area;            // Area of the element
    double TotalVolume;     // Total volume of mesh
    double TotalArea;       // Total area of mesh
    double Jacobian[3][3];   // Jacobian matrix
    double JacobianInv[3][3]; // Inverse of the Jacobian matrix
    double detJ;            // Determinant of the Jacobian
    double normal[3];       // Normal to the surface
} gelement3D;

// Function prototypes
void GetElement3D(const mesh* this, gelement3D* K, int i);
void Vertex_data(const mesh* this, vertex* Ver, int index);
void handleTriangle(const mesh* this, gelement3D* K, int index);
void handleQuadrangle(const mesh* this, gelement3D* K, int index);
void handleTetrahedron(const mesh* this, gelement3D* K, int index);
void handleHexahedron(const mesh* this, gelement3D* K, int index);
void handlePyramid(const mesh* this, gelement3D* K, int index);
void handlePrism(const mesh* this, gelement3D* K, int index);

double TetrahedronVolume(double *A, double *B, double *C, double *D);
int getBaseOffset(const mesh* this, int elementType);

void computeJacobianAndInverse(vertex Ver[], int Nvertex, int elementType, double J[3][3], double J_inv[3][3], double *det_J, double G[3][8], double xi, double eta, double zeta);
void computeSurfaceJacobianAndNormal(vertex Ver[], int nFaceNodes, double gradN_ref_xi[], double gradN_ref_eta[], double *J_surf, double normal[3], double tangents[2][3]);
#endif // GELEMENT3D_H