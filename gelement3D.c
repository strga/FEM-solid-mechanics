#include "gelement3D.h"
#include "mesh.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #define TRIAG_DATA
// #define QUAD_DATA
// #define TETRA_DATA
// #define HEX_DATA
// #define PYRAMID_DATA
// #define PRISM_DATA
// #define VOL_COMPUTE

// --------------------------------------------------------------------

void GetElement3D(const mesh* this, gelement3D* K, int i) {

    K->type = -1;

    int startQuad = this->NTriangles;
    int startTetra = startQuad + this->NQuads;
    int startHex = startTetra + this->NTetrahedra;
    int startPrism = startHex + this->NHexahedra;
    int startPyramid = startPrism + this->NWedges;

    if (i < startQuad) {
        handleTriangle(this, K, i);
        // printf("Element %d is a Triangle, K->type = %d\n", i, K->type);
        
    } else if (i >= startQuad && i < startTetra) {
        handleQuadrangle(this, K, i - startQuad);
        // printf("Element %d is a Quadrangle, K->type = %d\n", i, K->type);

    } else if (i >= startTetra && i < startHex) {
        handleTetrahedron(this, K, i - startTetra);
        // printf("Element %d is a Tetrahedron, K->type = %d\n", i, K->type);
        
    } else if (i >= startHex && i < startPrism) {
        handleHexahedron(this, K, i - startHex);
        // printf("Element %d is a Hexahedron, K->type = %d\n", i, K->type);
        
    } else if (i >= startPrism && i < startPyramid) {
        handlePrism(this, K, i - startPrism);
        // printf("Element %d is a Prism, K->type = %d\n", i, K->type);

    } else if (i >= startPyramid && i < this->NElements) {
        handlePyramid(this, K, i - startPyramid);
        // printf("Element %d is a Pyramid, K->type = %d\n", i, K->type);

    } else {
        
        fprintf(stderr, "Element index %d is out of bounds (Gelement3D)\n", i);
        exit(EXIT_FAILURE);
    }
}

int getBaseOffset(const mesh* this, int elementType) {
    int offset = 0;

    if (elementType > 0) offset += this->NTriangles * TRIANGLE_VERTEX;
    if (elementType > 1) offset += this->NQuads * QUADRANGLE_VERTEX;
    if (elementType > 2) offset += this->NTetrahedra * TETRAHEDRON_VERTEX;
    if (elementType > 3) offset += this->NHexahedra * HEXAHEDRON_VERTEX;
    if (elementType > 4) offset += this->NWedges * PRISM_VERTEX;
    if (elementType > 5) offset += this->NPyramids * PYRAMID_VERTEX;

    return offset;
}

void handleTriangle(const mesh* this, gelement3D* K, int index) {
    
    K->type = GMSH_TRIANGLE;
    K->Nvertex = TRIANGLE_VERTEX;
    // printf("handleTriangle called: K->type = %d, index = %d\n", K->type, index);
    double V[3], U[3], W[3];
    int baseOffset = getBaseOffset(this, 0);

    for (int k = 0; k < TRIANGLE_VERTEX; k++) {
        int vertex_index = this->ElementListOfVertices[baseOffset + TRIANGLE_VERTEX * index + k];
        K->idxNode[k] = vertex_index;
        K->Ver[k].x = this->x[vertex_index];
        K->Ver[k].y = this->y[vertex_index];
        K->Ver[k].z = this->z[vertex_index];
        // K->Ver[k].mark = this->ElementMark[vertex_index];

        #ifdef TRIAG_DATA
        if (k == 0) {
            printf("Triangle %d\n", index + 1);
        }
        printf("Vertex index: %d \t Index computed: %d\n", vertex_index + 1, 3 * index + k + 1);
        if (k == 2) {
            printf("\n");
        }
        #endif
    }

    // Calculate area using the cross product of vectors from the triangle vertices
    V[0] = K->Ver[1].x - K->Ver[0].x;
    V[1] = K->Ver[1].y - K->Ver[0].y;
    V[2] = K->Ver[1].z - K->Ver[0].z;

    U[0] = K->Ver[2].x - K->Ver[0].x;
    U[1] = K->Ver[2].y - K->Ver[0].y;
    U[2] = K->Ver[2].z - K->Ver[0].z;

    W[0] = V[1] * U[2] - V[2] * U[1];
    W[1] = V[2] * U[0] - V[0] * U[2];
    W[2] = V[0] * U[1] - V[1] * U[0];

    K->area = 0.5 * sqrt( W[0] * W[0] + W[1] * W[1] + W[2] * W[2] );
    K->TotalArea += K->area;
    // computeJacobianAndInverse(K->Ver, K->Nvertex, K->type, K->Jacobian, K->JacobianInv, &K->detJ);
    double invLength = 1.0 / (2.0 * K->area);
    K->normal[0] = W[0] * invLength;
    K->normal[1] = W[1] * invLength;
    K->normal[2] = W[2] * invLength;
}

void handleQuadrangle(const mesh* this, gelement3D* K, int index) {

    K->type = GMSH_QUADRANGLE;
    K->Nvertex = QUADRANGLE_VERTEX;
    // printf("handleQuadrangle called: K->type = %d, index = %d\n", K->type, index);
    double V1[3], V2[3], V3[3], V4[3], W1[3], W2[3];
    double area1, area2;
    int baseOffset = getBaseOffset(this, 1);

    for (int k = 0; k < QUADRANGLE_VERTEX; k++) {
        int vertex_index = this->ElementListOfVertices[baseOffset + QUADRANGLE_VERTEX * index + k];
        K->idxNode[k] = vertex_index;
        K->Ver[k].x = this->x[vertex_index];
        K->Ver[k].y = this->y[vertex_index];
        K->Ver[k].z = this->z[vertex_index];

        #ifdef QUAD_DATA
        if (k == 0) {
            printf("Quadrangle %d\n", index + 1);
        }
        printf("Vertex index: %d \t Index computed: %d\n", vertex_index + 1, 3 * quadIndex + 4 * index + k + 1);
        if (k == 2) {
            printf("\n");
        }
        #endif
    }

    // Vectors for two adjacent triangles (diagonal 1-3 split)
    V1[0] = K->Ver[1].x - K->Ver[0].x;
    V1[1] = K->Ver[1].y - K->Ver[0].y;
    V1[2] = K->Ver[1].z - K->Ver[0].z;

    V2[0] = K->Ver[3].x - K->Ver[0].x;
    V2[1] = K->Ver[3].y - K->Ver[0].y;
    V2[2] = K->Ver[3].z - K->Ver[0].z;

    V3[0] = K->Ver[2].x - K->Ver[1].x;
    V3[1] = K->Ver[2].y - K->Ver[1].y;
    V3[2] = K->Ver[2].z - K->Ver[1].z;

    V4[0] = K->Ver[3].x - K->Ver[1].x;
    V4[1] = K->Ver[3].y - K->Ver[1].y;
    V4[2] = K->Ver[3].z - K->Ver[1].z;

    // Cross products for two triangles
    W1[0] = V1[1] * V2[2] - V1[2] * V2[1];
    W1[1] = V1[2] * V2[0] - V1[0] * V2[2];
    W1[2] = V1[0] * V2[1] - V1[1] * V2[0];

    W2[0] = V3[1] * V4[2] - V3[2] * V4[1];
    W2[1] = V3[2] * V4[0] - V3[0] * V4[2];
    W2[2] = V3[0] * V4[1] - V3[1] * V4[0];

    // Area is the sum of two triangles
    area1 = 0.5 * sqrt(W1[0] * W1[0] + W1[1] * W1[1] + W1[2] * W1[2]);
    area2 = 0.5 * sqrt(W2[0] * W2[0] + W2[1] * W2[1] + W2[2] * W2[2]);

    K->area = area1 + area2;
    K->TotalArea += K->area;

    double W[3];
    W[0] = W1[0] + W2[0];
    W[1] = W1[1] + W2[1];
    W[2] = W1[2] + W2[2];

    double norm = sqrt(W[0] * W[0] + W[1] * W[1] + W[2] * W[2]); // Potential for improvement since quadrange can be nonplanar, then the result might nor be accurate since it is average of norms of 2 triangels

    if (norm < 1e-12) {
        fprintf(stderr, "Normal of a quadrangle is close to zero.\n");
        exit(EXIT_FAILURE);
    }

    K->normal[0] = W[0] / norm;
    K->normal[1] = W[1] / norm;
    K->normal[2] = W[2] / norm;

    // computeJacobianAndInverse(K->Ver, K->Nvertex, K->type, K->Jacobian, K->JacobianInv, &K->detJ);
}

void handleTetrahedron(const mesh* this, gelement3D* K, int index) {
    
    K->type = GMSH_TETRAHEDRON;
    K->Nvertex = TETRAHEDRON_VERTEX;
    // printf("handleTetrahedron called: K->type = %d, index = %d\n", K->type, index);
    double V[4][3];
    int baseOffset = getBaseOffset(this, 2);

    for (int k = 0; k < TETRAHEDRON_VERTEX; k++) {
        int vertex_index = this->ElementListOfVertices[baseOffset + TETRAHEDRON_VERTEX * index + k];
        K->idxNode[k] = vertex_index;
        K->Ver[k].x = this->x[vertex_index];
        K->Ver[k].y = this->y[vertex_index];
        K->Ver[k].z = this->z[vertex_index];
        // K->Ver[k].mark = this->ElementMark[vertex_index];

        V[k][0] = K->Ver[k].x;
        V[k][1] = K->Ver[k].y;
        V[k][2] = K->Ver[k].z;

        #ifdef TETRA_DATA
        if (k == 0) {
            printf("Tetrahedron %d\n", tetraIndex + index + 1);
        }
        printf("Vertex index: %d \t Index computed: %d\n", vertex_index + 1, 3 * tetraIndex + 4 * index + k );
        if (k == 3) {
            printf("\n");
        }
        #endif
    }

    K->volume = TetrahedronVolume(V[0], V[1], V[2], V[3]);
    K->TotalVolume += K->volume;
    // computeJacobianAndInverse(K->Ver, K->Nvertex, K->type, K->Jacobian, K->JacobianInv, &K->detJ);
    // printf("Tetrahedron volume: %g\n", K->volume);
}

void handleHexahedron(const mesh* this, gelement3D* K, int index) {
    
    K->type = GMSH_HEXAHEDRON;
    K->Nvertex = HEXAHEDRON_VERTEX;
    // printf("handleHexahedron called: K->type = %d, index = %d\n", K->type, index);
    double V[8][3];
    int baseOffset = getBaseOffset(this, 3);

    for (int k = 0; k < HEXAHEDRON_VERTEX; k++) {
        int vertex_index = this->ElementListOfVertices[baseOffset + HEXAHEDRON_VERTEX * index + k];
        K->idxNode[k] = vertex_index;
        K->Ver[k].x = this->x[vertex_index];
        K->Ver[k].y = this->y[vertex_index];
        K->Ver[k].z = this->z[vertex_index];
        // K->Ver[k].mark = this->ElementMark[vertex_index];

        V[k][0] = K->Ver[k].x;
        V[k][1] = K->Ver[k].y;
        V[k][2] = K->Ver[k].z;

        #ifdef HEX_DATA
        if (k == 0) {
            printf("Hexahedron %d\n", k + index + 1);
        }
        printf("Vertex index: %d \t Index computed: %d\n", vertex_index + 1, 8 * k + 4 * index + k );
        if (k == 8) {
            printf("\n");
        }
        #endif
    }

    // printf("Hexahedron %d: %d %d %d %d %d %d %d %d\n", index, K->idxNode[0], K->idxNode[1], K->idxNode[2], K->idxNode[3], K->idxNode[4], K->idxNode[5], K->idxNode[6], K->idxNode[7] );
    // for (int k = 0; k < HEXAHEDRON_VERTEX; k++)
    // {
    //     printf("Coords_idx %d: %g, %g, %g \n", K->idxNode[k], K->Ver[k].x, K->Ver[k].y, K->Ver[k].z);
    // }

    double volume = 0.0;
    volume += TetrahedronVolume(V[0], V[1], V[3], V[4]);
    volume += TetrahedronVolume(V[1], V[2], V[3], V[6]);
    volume += TetrahedronVolume(V[1], V[6], V[7], V[4]);
    volume += TetrahedronVolume(V[3], V[6], V[4], V[7]);
    volume += TetrahedronVolume(V[1], V[3], V[4], V[6]);

    K->volume = volume;
    K->TotalVolume += volume;
    // printf("Hexahedron before computeJacobianAndInverse\n");
    // printf("Hexahedron volume: %g\n", K->volume);
    // computeJacobianAndInverse(K->Ver, K->Nvertex, K->type, K->Jacobian, K->JacobianInv, &K->detJ);
}

void handlePyramid(const mesh* this, gelement3D* K, int index) {

    K->type = GMSH_PYRAMID;
    K->Nvertex = PYRAMID_VERTEX;
    // printf("handlePyramid called: K->type = %d, index = %d\n", K->type, index);
    double V[5][3];
    int baseOffset = getBaseOffset(this, 5);

    for (int k = 0; k < PYRAMID_VERTEX; k++) {
        int vertex_index = this->ElementListOfVertices[baseOffset + PYRAMID_VERTEX * index + k];
        K->idxNode[k] = vertex_index;
        K->Ver[k].x = this->x[vertex_index];
        K->Ver[k].y = this->y[vertex_index];
        K->Ver[k].z = this->z[vertex_index];
        // K->Ver[k].mark = this->ElementMark[vertex_index];

        V[k][0] = K->Ver[k].x;
        V[k][1] = K->Ver[k].y;
        V[k][2] = K->Ver[k].z;

        #ifdef PYRAMID_DATA
        if (k == 0) {
            printf("Pyramid %d\n", pyramidIndex + index + 1);
        }
        printf("Vertex index: %d \t Index computed: %d\n", vertex_index + 1, baseOffset + PYRAMID_VERTEX * index + k);
        if (k == 5) {
            printf("\n");
        }
        #endif
    }

    double volume = 0.0;
    volume += TetrahedronVolume(V[0], V[1], V[2], V[4]);
    volume += TetrahedronVolume(V[0], V[2], V[3], V[4]);
    
    K->volume = volume;
    K->TotalVolume += K->volume;
    // computeJacobianAndInverse(K->Ver, K->Nvertex, K->type, K->Jacobian, K->JacobianInv, &K->detJ);
    // printf("Pyramid volume: %g\n\n", K->volume);
}

void handlePrism(const mesh* this, gelement3D* K, int index) {

    K->type = GMSH_WEDGE;
    K->Nvertex = PRISM_VERTEX;
    // printf("handlePrism called: K->type = %d, index = %d\n", K->type, index);
    double V[6][3];
    int baseOffset = getBaseOffset(this, 4);

    for (int k = 0; k < PRISM_VERTEX; k++) {
        int vertex_index = this->ElementListOfVertices[baseOffset + PRISM_VERTEX * index + k];
        K->idxNode[k] = vertex_index;
        K->Ver[k].x = this->x[vertex_index];
        K->Ver[k].y = this->y[vertex_index];
        K->Ver[k].z = this->z[vertex_index];

        V[k][0] = K->Ver[k].x;
        V[k][1] = K->Ver[k].y;
        V[k][2] = K->Ver[k].z;

        #ifdef PRISM_DATA
        if (k == 0) {
            printf("Prism / Wedge %d\n", prismIndex + index + 1);
        }
        printf("Vertex index: %d \t Index computed: %d\n", vertex_index + 1, 4 * prismIndex + 6 * index + k + 1);
        if (k == 6) {
            printf("\n");
        }
        #endif
    }

    double volume = 0.0;
    volume += TetrahedronVolume(V[0], V[1], V[2], V[3]);
    volume += TetrahedronVolume(V[1], V[4], V[3], V[2]);
    volume += TetrahedronVolume(V[2], V[3], V[4], V[5]);

    K->volume = volume;
    K->TotalVolume += K->volume;
    // computeJacobianAndInverse(K->Ver, K->Nvertex, K->type, K->Jacobian, K->JacobianInv, &K->detJ);
    // printf("Prism volume: %g\n\n", K->volume);
}

double TetrahedronVolume(double *A, double *B, double *C, double *D) {
    double v1[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
    double v2[3] = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};
    double v3[3] = {D[0] - A[0], D[1] - A[1], D[2] - A[2]};

    double det = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
               - v1[1] * (v2[0] * v3[2] - v2[2] * v3[0])
               + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);

    double Vol_tetra = fabs(det) / 6.0;

    #ifdef VOL_COMPUTE
    printf("Vertices: A(%f, %f, %f), B(%f, %f, %f), C(%f, %f, %f), D(%f, %f, %f)\n", 
       A[0], A[1], A[2], B[0], B[1], B[2], C[0], C[1], C[2], D[0], D[1], D[2]);
    printf("v1 = (%f, %f, %f)\n", v1[0], v1[1], v1[2]);
    printf("v2 = (%f, %f, %f)\n", v2[0], v2[1], v2[2]);
    printf("v3 = (%f, %f, %f)\n", v3[0], v3[1], v3[2]);
    printf("Det: %f\n", det);
    printf("Sub_volume: %g\n", Vol_tetra);
    #endif

    return Vol_tetra;
}

void computeJacobianAndInverse(vertex Ver[], int Nvertex, int elementType, double J[3][3], double J_inv[3][3], double *det_J, double G[3][8], double xi, double eta, double zeta) {
    
    // static int k = 0;  // Persistent counter across function calls

    // printf("Execution count: k = %d\n", k);
    // printf("NEW ELEMENT:\n");
    // k++;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            J[i][j] = 0.0;
            J_inv[i][j] = 0.0;
        }
        for (int j = 0; j < Nvertex; j++) {
            G[i][j] = 0.0;
        }
    }
    
    double gradN_ref[3][8] = {0};
    GetReferenceGradients(elementType, gradN_ref[0], gradN_ref[1], gradN_ref[2], xi, eta, zeta);
    int nbBaseFn = GetNumBasisFunctions(elementType);
    // printf("nbBaseFn: %d\n", nbBaseFn);
    
    // for (int p = 0; p < nbBaseFn; p++)
    // {
    //     printf("Vertex %d: %g %g %g\n", p, Ver[p].x, Ver[p].y, Ver[p].z);
    // }
    // printf("\n");

    for (int i = 0; i < nbBaseFn; i++) {
        J[0][0] += gradN_ref[0][i] * Ver[i].x;
        J[1][0] += gradN_ref[0][i] * Ver[i].y;
        J[2][0] += gradN_ref[0][i] * Ver[i].z;

        // printf("Inside index = %d \t gradN_ref[0][index] = %f \t idx_i = %d ( %g %g %g )\n", i, gradN_ref[0][i], i, Ver[i].x, Ver[i].y, Ver[i].z);

        J[0][1] += gradN_ref[1][i] * Ver[i].x;
        J[1][1] += gradN_ref[1][i] * Ver[i].y;
        J[2][1] += gradN_ref[1][i] * Ver[i].z;

        // printf("Inside index = %d \t gradN_ref[1][index] = %f \t idx_i = %d ( %g %g %g )\n", i, gradN_ref[0][i], i, Ver[i].x, Ver[i].y, Ver[i].z);

        J[0][2] += gradN_ref[2][i] * Ver[i].x;
        J[1][2] += gradN_ref[2][i] * Ver[i].y;
        J[2][2] += gradN_ref[2][i] * Ver[i].z;

        // printf("Inside index = %d \t gradN_ref[2][index] = %f \t idx_i = %d ( %g %g %g )\n", i, gradN_ref[0][i], i, Ver[i].x, Ver[i].y, Ver[i].z);
        // printf("\n");

        // printf("x = %g y = %g z = %g\n", Ver[i].x, Ver[i].y, Ver[i].z);
        // printf("gradN_ref[%d][%d] = %g\n", 0, i, gradN_ref[0][i]);
        // printf("gradN_ref[%d][%d] = %g\n", 1, i, gradN_ref[1][i]);
        // printf("gradN_ref[%d][%d] = %g\n", 2, i, gradN_ref[2][i]);
        
    }
    
    // // Print Jacobian matrix J
    // printf("Jacobian matrix J:\n");
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         printf("%12.5f ", J[i][j]);
    //     }
    //     printf("\n");
    // }

    if (elementType == GMSH_TRIANGLE || elementType == GMSH_QUADRANGLE)
    {
        J[0][2] = J[1][2] = J[2][2] = 1.0;
    }
    
    *det_J = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
           - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
           + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    
    // printf("\n");
    // printf("det(J) = %12.5f\n", *det_J);

    if (fabs(*det_J) < 1e-12) {
        fprintf(stderr, "Jacobian determinant is near zero. Element may be degenerate.\n");
        exit(EXIT_FAILURE);
    }

    double inv_det_J = 1.0 / (*det_J);

    J_inv[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv_det_J;
    J_inv[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * inv_det_J;
    J_inv[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv_det_J;

    J_inv[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * inv_det_J;
    J_inv[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv_det_J;
    J_inv[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * inv_det_J;

    J_inv[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv_det_J;
    J_inv[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * inv_det_J;
    J_inv[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv_det_J;

    // // Print inverse Jacobian matrix J_inv
    // printf("\n");
    // printf("Inverse Jacobian matrix J_inv:\n");
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         printf("%12.5f ", J_inv[i][j]);
    //     }
    //     printf("\n");
    // }

    for (int i = 0; i < nbBaseFn; i++) {
        for (int k = 0; k < 3; k++) { // physical coordinate (x, y, z)
            G[k][i] = 0.0;
            for (int l = 0; l < 3; l++) { // reference coordinate (ξ, η, ζ)
                G[k][i] += J_inv[l][k] * gradN_ref[l][i];
            }
        }
    }
 
    // Print transformed gradients in physical space (G matrix)
    // printf("\n");
    // printf("Gradients in Physical Space (G matrix):\n");
    // for (int i = 0; i < nbBaseFn; i++) {
    //     printf("G[%d] = (%12.5f, %12.5f, %12.5f)\n",
    //            i, G[0][i], G[1][i], G[2][i]);
    // }
}

void computeSurfaceJacobianAndNormal(vertex Ver[], int nFaceNodes, double gradN_ref_xi[], double gradN_ref_eta[], double *J_surf, double normal[3], double tangents[2][3]) {
    double PxH[3] = {0.0}, PyH[3] = {0.0};

    // Local frame on face
    for (int i = 0; i < nFaceNodes; i++) {
        PxH[0] += gradN_ref_xi[i] * Ver[i].x;
        PxH[1] += gradN_ref_xi[i] * Ver[i].y;
        PxH[2] += gradN_ref_xi[i] * Ver[i].z;

        PyH[0] += gradN_ref_eta[i] * Ver[i].x;
        PyH[1] += gradN_ref_eta[i] * Ver[i].y;
        PyH[2] += gradN_ref_eta[i] * Ver[i].z;
    }

    // Normal = Px × Py
    double n[3];
    n[0] = PxH[1] * PyH[2] - PxH[2] * PyH[1];
    n[1] = PxH[2] * PyH[0] - PxH[0] * PyH[2];
    n[2] = PxH[0] * PyH[1] - PxH[1] * PyH[0];

    double norm_n = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    if (norm_n < 1e-12) {
        fprintf(stderr, "Degenerate face detected (surface Jacobian near zero)\n");
        exit(EXIT_FAILURE);
    }

    *J_surf = norm_n;
    normal[0] = n[0] / norm_n;
    normal[1] = n[1] / norm_n;
    normal[2] = n[2] / norm_n;

    // Optional: tangents on surface
    if (tangents) {
        tangents[0][0] = PxH[0];
        tangents[0][1] = PxH[1];
        tangents[0][2] = PxH[2];

        tangents[1][0] = PyH[0];
        tangents[1][1] = PyH[1];
        tangents[1][2] = PyH[2];
    }
}



