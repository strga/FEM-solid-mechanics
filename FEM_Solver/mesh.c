#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "mesh.h"
#include "wlog.h"
#include "gmshtools.h"
#include "integral.h"
#include "triplet.h"
#include "Triplet_1.h"
#include "gelement.h"
#include "gelement3D.h"
#include "Mesh_GetBndrSide.h"
#include "saving.h"

// #define RHS_INTEGRAL
// #define ROBIN_VALUES
// #define LOAD_VECTOR
// #define Bmatrix
// #define Dmatrix

/* ------------------------------------------------------------------------------------------------- */
void Mesh_Allocate(mesh *this)
{
    int num;
    debuglog("Debug[MESH]:\t\t Allocating mesh.\n");

    this->x                   = malloc(4 * this->NPoints * sizeof(double)); MemCheck(this->x, __func__);
    this->y                   = this->x + this->NPoints;

    if (this->dim == 3) {
        this->z               = this->y + this->NPoints;
    }
    else {
        this->z               = NULL;
    }
    
    this->NodeDegree          = malloc(this->NPoints * sizeof(short));  MemCheck(this->NodeDegree,  __func__);
    this->PointMark           = malloc(this->NPoints * sizeof(short));  MemCheck(this->PointMark, __func__);
    
    /* Allocating ELEMENTS  */
    debuglog("Debug[ELEMENTS]:\t Allocating memory (T%d, Q%d).\n", this->NTriangles, this->NQuads);
    
    
    // num = 3 * this->NTriangles + 4 * this->NQuads;
    
    if (this->dim == 3)
    {
        num = 3 * this->NTriangles   // 2D triangles (for surfaces or interfaces)
        + 4 * this->NQuads           // 2D quadrangles (for surfaces or interfaces)
        + 4 * this->NTetrahedra      // 3D tetrahedra
        + 8 * this->NHexahedra       // 3D hexahedra
        + 6 * this->NWedges          // 3D wedges (prisms)
        + 5 * this->NPyramids;       // 3D pyramids
    }
    else
    {
        num = 3 * this->NTriangles + 4 * this->NQuads;
    }
    
    this->ElementListOfVertices = malloc(num * sizeof(int));  MemCheck(this->ElementListOfVertices, __func__);
    this->ElementMark           = malloc((this->NElements) * sizeof(short)); MemCheck(this->ElementMark, __func__);
    debuglog("Debug[ELEMENTS]:\t Allocating memory. \t\t\t (done)\n");
    
    /* Allocating Edges Edges_Allocate(&this->sideh); */
    //int     *EdgeVerA, *EdgeVerB;

    this->EdgeVerA   = malloc(2 * this->NBEdges * sizeof(int));          MemCheck(this->EdgeVerA, __func__);
    this->EdgeVerB   = this->EdgeVerA + this->NBEdges;
    this->EdgeMark   = malloc(this->NBEdges * sizeof(short));        MemCheck(this->EdgeMark, __func__);
    
    //this->Mesh_FillPatches = malloc((this->NElements) * sizeof(int));

    debuglog("Debug[MESH]:\t\t Allocating mesh. \t\t(done)\n");
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Free(mesh *this)
{
    //int ii;
    debuglog("Debug[MESH]:\t\t Deallocating mesh. \n");
    
    free(this->x); 
    free(this->NodeDegree);
    free(this->PointMark);
    this->PointMark = NULL;

    free(this->ElementListOfVertices);
    free(this->ElementMark);
    this->ElementListOfVertices = NULL;
    this->ElementMark           = NULL;

    free(this->EdgeVerA);
    free(this->EdgeMark);


    //free(this->PatchListOfElements);
    //free(this->PatchW);
    //free(this->PatchListOfElements);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Load(mesh *this, const char *fname)
{
    Mesh_ImportMSH(this, fname); 
    PointMark(this);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_ImportElementsMSH(mesh *this, const char *fname);
/* ------------------------------------------------------------------------------------------------- */
void Mesh_ImportMSH(mesh *this, const char *fname)
{
    FILE *fid;
    //int   np; 
    //int   j, i;
    gmshfile grid;

    // Initialize mesh element counts
    this->NTriangles = this->NQuads =  0;
    this->NTetrahedra = this->NHexahedra = this->NWedges = this->NPyramids = 0;

    gmshfile_read(&grid, fname);
    int dim = (grid.NTetrahedra > 0 || grid.NHexahedra > 0 || grid.NWedges > 0 || grid.NPyramids > 0) ? 3 : 2; // Determine dimension
    this->dim = dim;

    // gmshfile_read(&grid, fname);

    this->NPoints           = grid.NPoints;
    // this->NElements         = grid.NTriangles + grid.NQuads;
    this->NElements         = grid.NTriangles + grid.NQuads + grid.NTetrahedra + grid.NHexahedra + grid.NWedges + grid.NPyramids;
    this->NTriangles        = grid.NTriangles;
    this->NQuads            = grid.NQuads;
    this->NTetrahedra       = grid.NTetrahedra;
    this->NHexahedra        = grid.NHexahedra;
    this->NWedges           = grid.NWedges;
    this->NPyramids         = grid.NPyramids;

    this->NBEdges           = grid.NBndrEdges;

    this->NTetrahedra   = grid.NTetrahedra;
    this->NHexahedra    = grid.NHexahedra;
    this->NWedges       = grid.NWedges;
    this->NPyramids     = grid.NPyramids;

    Mesh_Allocate(this);

    fid = fopen(fname, "r"); FileCheck(fid, "Error opening MSH file %s", fname);
    gmsh_getNodesCoordinates(fid, this->x, this->y, this->z, dim);
    fclose(fid);
    Mesh_ImportElementsMSH(this, fname);

    writelog("Init[GMESH/MSH]:\t Load mesh from %s  \t\t\t\t\t\t\t [ok]\n", fname);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_ImportElementsMSH(mesh *this, const char *fname)
{
    FILE *fid;
    char buf[BUFFER_SIZE + 1];
    int idxQuad, idxTriangle, idxEdge, idxTetra, idxHex, idxWedge, idxPyramid;
    int startQuads, startTetra, startHex, startWedge, startPyramid;
    int nall, i;// idx, j;
    gmshline current;
    //int *pomWritten;
    //int SumNodeDegree;

    fid = fopen(fname, "r"); FileCheck(fid, "Elements::Load file %s\n", fname);

    idxQuad = idxTriangle = idxEdge = idxTetra = idxHex = idxWedge = idxPyramid = 0;
    
    // startQuads = this->NTriangles * 3;
    // startTetra = startQuads + this->NQuads * 4;
    // startHex = startTetra + this->NTetrahedra * 4;
    // startPyramid = startHex + this->NWedges * 6;
    // startWedge = startPyramid + this->NHexahedra * 8;

    startQuads = this->NTriangles * 3;
    startTetra = startQuads + this->NQuads * 4;
    startHex = startTetra + this->NTetrahedra * 4;
    startWedge = startHex + this->NHexahedra * 8;
    startPyramid = startWedge + this->NWedges * 6;

    // printf("Offsets: Quads: %d, Tetra: %d, Hex: %d, Wedge: %d, Pyramid: %d\n", startQuads, startTetra, startHex, startWedge, startPyramid);
    // getchar();

    do fgets(buf, BUFFER_SIZE, fid); while (strncasecmp(buf, "$Elements", 9) != 0 && feof(fid) == 0 );
    fgets(buf, BUFFER_SIZE, fid); 
    sscanf(buf, "%d", &nall);
    for(i = 0; i < nall; i++)
    {
        fgets(buf, BUFFER_SIZE, fid); 
        gmshline_read(&current, buf);
        switch (current.etyp)
        {
        case GMSH_SEGMENT:
            this->EdgeMark[idxEdge]      = current.tagPhysical;
            this->EdgeVerA[idxEdge]      = current.ilist[0];
            this->EdgeVerB[idxEdge]      = current.ilist[1];
            idxEdge++;
            break;

        case GMSH_TRIANGLE:
            this->ElementListOfVertices[3 * idxTriangle + 0] = current.ilist[0];
            this->ElementListOfVertices[3 * idxTriangle + 1] = current.ilist[1];
            this->ElementListOfVertices[3 * idxTriangle + 2] = current.ilist[2];
            this->ElementMark[idxTriangle]      = current.tagPhysical;
            
               this->NodeDegree[current.ilist[0]]++;
               this->NodeDegree[current.ilist[1]]++;
               this->NodeDegree[current.ilist[2]]++;

            idxTriangle++;
            break;

        case GMSH_QUADRANGLE:
            this->ElementListOfVertices[startQuads + 4 * idxQuad + 0] = current.ilist[0];
            this->ElementListOfVertices[startQuads + 4 * idxQuad + 1] = current.ilist[1];
            this->ElementListOfVertices[startQuads + 4 * idxQuad + 2] = current.ilist[2];
            this->ElementListOfVertices[startQuads + 4 * idxQuad + 3] = current.ilist[3];
            this->ElementMark[this->NTriangles + idxQuad]           = current.tagPhysical;

               this->NodeDegree[current.ilist[0]]++;
               this->NodeDegree[current.ilist[1]]++;
               this->NodeDegree[current.ilist[2]]++;
               this->NodeDegree[current.ilist[3]]++;
            
            idxQuad++;
            break;

        case GMSH_TETRAHEDRON:
            this->ElementListOfVertices[startTetra + 4 * idxTetra + 0] = current.ilist[0];
            this->ElementListOfVertices[startTetra + 4 * idxTetra + 1] = current.ilist[1];
            this->ElementListOfVertices[startTetra + 4 * idxTetra + 2] = current.ilist[2];
            this->ElementListOfVertices[startTetra + 4 * idxTetra + 3] = current.ilist[3];
            this->ElementMark[this->NTriangles + this->NQuads + idxTetra] = current.tagPhysical;

                this->NodeDegree[current.ilist[0]]++;
                this->NodeDegree[current.ilist[1]]++;
                this->NodeDegree[current.ilist[2]]++;
                this->NodeDegree[current.ilist[3]]++;

            // printf("Vertex1: %d \t ilist_0: %d\n",startTetra + 4 * idxTetra + 0, current.ilist[0] + 1);
            // printf("Vertex2: %d \t ilist_1: %d\n",startTetra + 4 * idxTetra + 1, current.ilist[1] + 1);
            // printf("Vertex3: %d \t ilist_2: %d\n",startTetra + 4 * idxTetra + 2, current.ilist[2] + 1);
            // printf("Vertex4: %d \t ilist_3: %d\n",startTetra + 4 * idxTetra + 3, current.ilist[3] + 1); 
            // printf("Tetra_index: %d\n", idxTetra);
            // printf("\n");

            // ALL OK

            idxTetra++;
            break;

        case GMSH_HEXAHEDRON:
            this->ElementListOfVertices[startHex + 8 * idxHex + 0] = current.ilist[0];
            this->ElementListOfVertices[startHex + 8 * idxHex + 1] = current.ilist[1];
            this->ElementListOfVertices[startHex + 8 * idxHex + 2] = current.ilist[2];
            this->ElementListOfVertices[startHex + 8 * idxHex + 3] = current.ilist[3];
            this->ElementListOfVertices[startHex + 8 * idxHex + 4] = current.ilist[4];
            this->ElementListOfVertices[startHex + 8 * idxHex + 5] = current.ilist[5];
            this->ElementListOfVertices[startHex + 8 * idxHex + 6] = current.ilist[6];
            this->ElementListOfVertices[startHex + 8 * idxHex + 7] = current.ilist[7];
            this->ElementMark[this->NTriangles + this->NQuads + this->NTetrahedra + idxHex] = current.tagPhysical;

                // Increment NodeDegree for each vertex involved in the Hexahedron
                this->NodeDegree[current.ilist[0]]++;
                this->NodeDegree[current.ilist[1]]++;
                this->NodeDegree[current.ilist[2]]++;
                this->NodeDegree[current.ilist[3]]++;
                this->NodeDegree[current.ilist[4]]++;
                this->NodeDegree[current.ilist[5]]++;
                this->NodeDegree[current.ilist[6]]++;
                this->NodeDegree[current.ilist[7]]++;

            idxHex++;
            break;

        case GMSH_WEDGE:
            this->ElementListOfVertices[startWedge + 6 * idxWedge + 0] = current.ilist[0];
            this->ElementListOfVertices[startWedge + 6 * idxWedge + 1] = current.ilist[1];
            this->ElementListOfVertices[startWedge + 6 * idxWedge + 2] = current.ilist[2];
            this->ElementListOfVertices[startWedge + 6 * idxWedge + 3] = current.ilist[3];
            this->ElementListOfVertices[startWedge + 6 * idxWedge + 4] = current.ilist[4];
            this->ElementListOfVertices[startWedge + 6 * idxWedge + 5] = current.ilist[5];
            this->ElementMark[this->NTriangles + this->NQuads + this->NTetrahedra + this->NHexahedra + idxWedge] = current.tagPhysical;

                // Increment NodeDegree for each vertex involved in the Wedge
                this->NodeDegree[current.ilist[0]]++;
                this->NodeDegree[current.ilist[1]]++;
                this->NodeDegree[current.ilist[2]]++;
                this->NodeDegree[current.ilist[3]]++;
                this->NodeDegree[current.ilist[4]]++;
                this->NodeDegree[current.ilist[5]]++;

            idxWedge++;
            break;

        case GMSH_PYRAMID:
            this->ElementListOfVertices[startPyramid + 5 * idxPyramid + 0] = current.ilist[0];
            this->ElementListOfVertices[startPyramid + 5 * idxPyramid + 1] = current.ilist[1];
            this->ElementListOfVertices[startPyramid + 5 * idxPyramid + 2] = current.ilist[2];
            this->ElementListOfVertices[startPyramid + 5 * idxPyramid + 3] = current.ilist[3];
            this->ElementListOfVertices[startPyramid + 5 * idxPyramid + 4] = current.ilist[4];
            this->ElementMark[this->NTriangles + this->NQuads + this->NTetrahedra + this->NHexahedra + this->NWedges + idxPyramid] = current.tagPhysical;

                // Incrementing the NodeDegree for each vertex involved in the Pyramid
                this->NodeDegree[current.ilist[0]]++;
                this->NodeDegree[current.ilist[1]]++;
                this->NodeDegree[current.ilist[2]]++;
                this->NodeDegree[current.ilist[3]]++;
                this->NodeDegree[current.ilist[4]]++;

                // printf("Vertex1: %d \t ilist_0: %d\n",startPyramid + 5 * idxPyramid + 0, current.ilist[0] + 1);
                // printf("Vertex2: %d \t ilist_1: %d\n",startPyramid + 5 * idxPyramid + 1, current.ilist[1] + 1);
                // printf("Vertex3: %d \t ilist_2: %d\n",startPyramid + 5 * idxPyramid + 2, current.ilist[2] + 1);
                // printf("Vertex4: %d \t ilist_3: %d\n",startPyramid + 5 * idxPyramid + 3, current.ilist[3] + 1);
                // printf("Vertex5: %d \t ilist_4: %d\n",startPyramid + 5 * idxPyramid + 4, current.ilist[4] + 1);
                // printf("Pyramid_index: %d\n", idxTetra);
                // printf("\n");

            idxPyramid++;
            break;

        default:
            writelog("ERROR:Wrong element!\n");
            break;
        }
    }

    fgets(buf, BUFFER_SIZE, fid);
    if (strncasecmp(buf, "$EndElements", 12)!=0)
        CriticalErrorMsg("Mesh:ReadElementsFromFileMSH\tIncorrect mesh format (end elements)\n");
    fclose(fid);

}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_FillPatches(mesh *this)
{
    
    int i, j, idx;
    int *PatchPom;
    int SumNodeDegree;

    // Calculate the sum of node degrees
    SumNodeDegree = 0;
    for (i = 0; i < this->NPoints; i++) {
        SumNodeDegree += this->NodeDegree[i];
    }

    // Allocate memory for PatchW and initialize
    this->PatchW = malloc((this->NPoints + 1) * sizeof(int));
    if (this->PatchW == NULL) {
        // Handle the memory allocation error
        fprintf(stderr, "Error allocating memory for PatchW\n");
        return;
    }
    this->PatchW[0] = 0;
    for (i = 1; i < this->NPoints + 1; i++) {
        this->PatchW[i] = this->PatchW[i-1] + this->NodeDegree[i-1];
    }

    // Allocate memory for PatchListOfElements
    this->PatchListOfElements = malloc(SumNodeDegree * sizeof(int));
    if (this->PatchListOfElements == NULL) {
        // Handle the memory allocation error and free already allocated memory
        fprintf(stderr, "Error allocating memory for PatchListOfElements\n");
        free(this->PatchW);
        return;
    }
    PatchPom = malloc(this->NPoints * sizeof(int));
    if (PatchPom == NULL) {
        // Handle the memory allocation error and free already allocated memory
        fprintf(stderr, "Error allocating memory for PatchPom\n");
        free(this->PatchW);
        free(this->PatchListOfElements);
        return;
    }

    // Initialize PatchPom array
    for (i = 0; i < this->NPoints; i++) {
        PatchPom[i] = 0; 
    }

    // Fill PatchListOfElements
    for (i = 0; i < this->NElements; i++) {
        for (j = 0; j < 3; j++) { // Assuming 3 vertices per element for triangles
            idx = this->ElementListOfVertices[3 * i + j];
            this->PatchListOfElements[this->PatchW[idx] + PatchPom[idx]] = i;
            PatchPom[idx]++;
        }
    }

    // Free allocated memory for PatchPom
    free(PatchPom);
}
/* ------------------------------------------------------------------------------------------------- */
void PointMark(mesh *this)
{
    
    int i, iA, iB;

    for ( i = 0; i < this->NPoints; i++)
    {
        this->PointMark[i] = 0;
    }
    
    for ( i = 0; i < this->NBEdges; i++)
    {
        iA = this->EdgeVerA[i];
        iB = this->EdgeVerB[i];

        if (this->EdgeMark[i] != 0)
        {
            if (this->PointMark[iA] == 0)
            {
                this->PointMark[iA] = this->EdgeMark[i];
            }

            else
                this->PointMark[iA] = min(this->PointMark[iA], this->EdgeMark[i]);

        
            if (this->PointMark[iB] == 0)
            {
                this->PointMark[iB] = this->EdgeMark[i];
            }

            else
                this->PointMark[iB] = min(this->PointMark[iB], this->EdgeMark[i]);

        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void IdentifyEdgePoints(mesh *this, int **edgePoints)
{

    *edgePoints = (int *)malloc(this->NPoints * sizeof(int));

    if (*edgePoints == NULL)
    {
        fprintf(stderr, "Error allocating memory for edge points\n");
        return;
    }

    for (int i = 0; i < this->NPoints; i++)
    {

        (*edgePoints)[i] = (this->PointMark[i] != 0) ? 1 : 0; // if condition true then 1 elese 0
    }
}
/* ------------------------------------------------------------------------------------------------- */
void IdentifyPointMark1(mesh *this, int **isDirichlet)
{
    int *edgePoints;
    IdentifyEdgePoints(this, &edgePoints);

    *isDirichlet = (int *)malloc(this->NPoints * sizeof(int));

    if (*isDirichlet == NULL)
    {
        fprintf(stderr, "Error allocating memory for logical vector\n");

        return;
    }

    for (int i = 0; i < this->NPoints; i++)
    {
        if (edgePoints[i] && (this->PointMark[i] == 1))
        {
            (*isDirichlet)[i] = 1;
        }
        else
        {
            (*isDirichlet)[i] = 0;
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_FillPatches3D(mesh *this) {
    int i, j, idx;
    int *PatchPom;
    int SumNodeDegree = 0;

    // Calculate the sum of node degrees (total number of element-node connections)
    for (i = 0; i < this->NPoints; i++) {
        SumNodeDegree += this->NodeDegree[i];
    }

    // Allocate memory for PatchW
    this->PatchW = malloc((this->NPoints + 1) * sizeof(int));
    if (this->PatchW == NULL) {
        fprintf(stderr, "Error allocating memory for PatchW\n");
        return;
    }
    this->PatchW[0] = 0;
    for (i = 1; i <= this->NPoints; i++) {
        this->PatchW[i] = this->PatchW[i - 1] + this->NodeDegree[i - 1];
    }

    // Allocate memory for PatchListOfElements
    this->PatchListOfElements = malloc(SumNodeDegree * sizeof(int));
    if (this->PatchListOfElements == NULL) {
        fprintf(stderr, "Error allocating memory for PatchListOfElements\n");
        free(this->PatchW);
        return;
    }
    PatchPom = calloc(this->NPoints, sizeof(int));
    if (PatchPom == NULL) {
        fprintf(stderr, "Error allocating memory for PatchPom\n");
        free(this->PatchW);
        free(this->PatchListOfElements);
        return;
    }

    // Fill PatchListOfElements
    for (i = 0; i < this->NElements; i++) {
        // Get the number of vertices for the current element
        int nVertices = 0;
        switch (this->ElementMark[i]) {
            case GMSH_TETRAHEDRON: nVertices = 4; break;
            case GMSH_HEXAHEDRON:  nVertices = 8; break;
            case GMSH_WEDGE:       nVertices = 6; break;
            case GMSH_PYRAMID:     nVertices = 5; break;
            default:
                fprintf(stderr, "Unknown element type in Mesh_FillPatches3D\n");
                continue;
        }

        // Associate each vertex of the element with the element itself
        for (j = 0; j < nVertices; j++) {
            idx = this->ElementListOfVertices[nVertices * i + j]; // Node index
            this->PatchListOfElements[this->PatchW[idx] + PatchPom[idx]] = i; // Store the element
            PatchPom[idx]++;
        }
    }

    // Free allocated memory for PatchPom
    free(PatchPom);
}
/* ------------------------------------------------------------------------------------------------- */
void IdentifyFacePoints(mesh *this, int **facePoints) {
    *facePoints = (int *)malloc(this->NPoints * sizeof(int));
    if (*facePoints == NULL) {
        fprintf(stderr, "Error allocating memory for face points\n");
        return;
    }

    // A point is on the boundary if its PointMark is nonzero
    for (int i = 0; i < this->NPoints; i++) {
        (*facePoints)[i] = (this->PointMark[i] != 0) ? 1 : 0;
    }
}
/* ------------------------------------------------------------------------------------------------- */
void IdentifyPointMark3D(mesh *this, int **isDirichlet) {
    int *facePoints;
    IdentifyFacePoints(this, &facePoints);

    *isDirichlet = (int *)malloc(this->NPoints * sizeof(int));
    if (*isDirichlet == NULL) {
        fprintf(stderr, "Error allocating memory for Dirichlet points\n");
        free(facePoints);
        return;
    }

    for (int i = 0; i < this->NPoints; i++) {
        if (facePoints[i] && (this->PointMark[i] == DIRICHLET_BC)) {
            (*isDirichlet)[i] = 1;
        } else {
            (*isDirichlet)[i] = 0;
        }
    }

    free(facePoints);
}
/* ------------------------------------------------------------------------------------------------- */
void IdentifyFaceMark3D(mesh *this, int **faceMark) {
    *faceMark = (int *)malloc(this->NPoints * sizeof(int));

    for (int i = 0; i < this->NPoints; i++) {
        (*faceMark)[i] = (this->PointMark[i] != 0) ? this->PointMark[i] : 0;
    }
}
/* ------------------------------------------------------------------------------------------------- */
double CalculateIntegral(mesh *this, gelement K, double (*f)(double, double), IntegrationType type, int i)
{
    double integralValue = 0.0;
    double omega = 1.0 / 3.0; // Weight

    switch (type)
    {
        case INTEGRATION_VERTICES:
        {
            double ff[3] = { f(K.A[0], K.A[1]), f(K.B[0], K.B[1]), f(K.C[0], K.C[1]) };
            integralValue = K.dV * omega * ff[i];
            break;
        }
        
        case INTEGRATION_MIDPOINTS:
        {
            // Midpoints of each side of the triangle
            double midAB[2] = {(K.A[0] + K.B[0]) / 2.0, (K.A[1] + K.B[1]) / 2.0};
            double midBC[2] = {(K.B[0] + K.C[0]) / 2.0, (K.B[1] + K.C[1]) / 2.0};
            double midCA[2] = {(K.C[0] + K.A[0]) / 2.0, (K.C[1] + K.A[1]) / 2.0};

            // Function values at midpoints
            double fMidAB = f(midAB[0], midAB[1]);
            double fMidBC = f(midBC[0], midBC[1]);
            double fMidCA = f(midCA[0], midCA[1]);

            // Adjusted integral calculation for each midpoint
            double fMid[3] = {fMidAB, fMidBC, fMidCA};
            integralValue = K.dV * omega * fMid[i];
            break;
        }

        case INTEGRATION_TEZISTE:
        {
            // Centroid of the triangle
            double centroid[2] = {(K.A[0] + K.B[0] + K.C[0]) / 3.0, (K.A[1] + K.B[1] + K.C[1]) / 3.0};

            // Function value at the centroid
            double fCentroid = f(centroid[0], centroid[1]);

            // Integral value at the centroid for each vertex
            integralValue = K.dV * omega * fCentroid;
            break;
        }

        case INTEGRATION_GAUSSIAN:
        {
            // Centroid of the triangle, as the Gaussian integration point
            double centroidX = (K.A[0] + K.B[0] + K.C[0]) / 3.0;
            double centroidY = (K.A[1] + K.B[1] + K.C[1]) / 3.0;

            // The weight for 1-point Gaussian Quadrature in the area coordinate system for a triangle
            double weight = K.dV * 0.5; // Area of triangle times the weight of 1/2

            integralValue = weight * f(centroidX, centroidY);
            break;
        }

        default:
            fprintf(stderr, "Error: Unknown integration type.\n");
    }

    return integralValue;
}
/* ------------------------------------------------------------------------------------------------- */
void computeGradient(const gelement *K, double G[2][3])
{
  
    // Gradients of base funciton for triangular (reference) element
    double dPhi_xy[3][2] = {{-1, -1}, {1, 0}, {0, 1}};

    // Gradient matrix
    for (int i = 0; i < 2; ++i) //Loop over columns of G  (dx, dy)
    {
        for (int j = 0; j < 3; ++j) // Loop over rows of G
        {
            // Compute the transpose of gradient of elements ... dx + dy
            G[i][j] = dPhi_xy[j][0] * K->invBk[0][i] + dPhi_xy[j][1] * K->invBk[1][i];
        }
    }

    // printf("Inverse Jacobian (Bk):\n");
    // for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 2; j++) {
    //         printf("%f ", K->invBk[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("Computed Gradients (G):\n");
    // for (int i = 0; i < 2; i++) { // dx, dy
    //     for (int j = 0; j < 3; j++) { // Node gradients
    //         printf("G[%d][%d] = %f\n", i, j, G[i][j]);
    //     }
    // }
    // printf("\n");
}
/* ------------------------------------------------------------------------------------------------- */
void ComputeBasisFunctions(int elementType, double *N, double xi, double eta, double zeta, int nNodes) {
    switch (elementType) {
        case GMSH_TRIANGLE: // Triangle (3 Nodes, 2D)
            N[0] = 1 - xi - eta;
            N[1] = xi;
            N[2] = eta;

            break;

        case GMSH_QUADRANGLE: // Quadrilateral (4 Nodes, 2D) - sometimes N[2] and N[3] are swapped
            N[0] = (1 - xi) * (1 - eta);
            N[1] = xi * (1 - eta);
            N[2] = xi * eta;
            N[3] = (1 - xi) * eta;

            break;

         case GMSH_TETRAHEDRON: // Tetrahedron (4 Nodes, 3D)
            N[0] = 1 - xi - eta - zeta;
            N[1] = xi;
            N[2] = eta;
            N[3] = zeta;

            break;

        case GMSH_HEXAHEDRON: // Hexahedron (8 Nodes, 3D)
            N[0] = (1 - xi) * (1 - eta) * (1 - zeta);
            N[1] = xi * (1 - eta) * (1 - zeta);
            N[2] = xi * eta * (1 - zeta);
            N[3] = (1 - xi) * eta * (1 - zeta);
            N[4] = (1 - xi) * (1 - eta) * zeta;
            N[5] = xi * (1 - eta) * zeta;
            N[6] = xi * eta * zeta;
            N[7] = (1 - xi) * eta * zeta;

            break;

        case GMSH_WEDGE: // Prism (6 Nodes, 3D)
            N[0] = (1 - xi - eta) * (1 - zeta);
            N[1] = xi * (1 - zeta);
            N[2] = eta * (1 - zeta);
            N[3] = (1 - xi - eta) * zeta;
            N[4] = xi * zeta;
            N[5] = eta * zeta;

            break;

        case GMSH_PYRAMID: // Pyramid (5 Nodes, 3D)
            N[0] = ( 1 - xi - eta - zeta + ( xi * eta ) / ( 1 - zeta ) ) / 4.0;
            N[1] = ( 1 + xi - eta - zeta - ( xi * eta ) / ( 1 - zeta ) ) / 4.0;
            N[2] = ( 1 + xi + eta - zeta + ( xi * eta ) / ( 1 - zeta ) ) / 4.0;
            N[3] = ( 1 - xi + eta - zeta - ( xi * eta ) / ( 1 - zeta ) ) / 4.0;
            N[4] = zeta; // Apex

            break;

        default:
            fprintf(stderr, "Unknown element type: %d\n", elementType);
            exit(EXIT_FAILURE);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void GetReferenceGradients(int elementType, double *dN_xi, double *dN_eta, double *dN_zeta, double xi, double eta, double zeta) {
    switch (elementType) {
        case GMSH_TRIANGLE:
            dN_xi[0] = -1.0;   dN_eta[0] = -1.0;
            dN_xi[1] = 1.0;    dN_eta[1] = 0.0;
            dN_xi[2] = 0.0;    dN_eta[2] = 1.0;
            break;

        case GMSH_QUADRANGLE:
            dN_xi[0] = eta - 1;    dN_eta[0] = xi - 1;
            dN_xi[1] = 1 - eta;    dN_eta[1] = -xi;
            dN_xi[2] = eta;        dN_eta[2] = xi;
            dN_xi[3] = -eta;       dN_eta[3] = 1 - xi;
            break;

        case GMSH_TETRAHEDRON:
            dN_xi[0] = -1.0;   dN_eta[0] = -1.0;  dN_zeta[0] = -1.0;
            dN_xi[1] = 1.0;    dN_eta[1] = 0.0;   dN_zeta[1] = 0.0;
            dN_xi[2] = 0.0;    dN_eta[2] = 1.0;   dN_zeta[2] = 0.0;
            dN_xi[3] = 0.0;    dN_eta[3] = 0.0;   dN_zeta[3] = 1.0;
            break;

            case GMSH_HEXAHEDRON:
            dN_xi[0] = - (1 - eta) * (1 - zeta);    dN_eta[0] = - (1 - xi) * (1 - zeta);    dN_zeta[0] = - (1 - xi) * (1 - eta);
            dN_xi[1] = (1 - eta) * (1 - zeta);      dN_eta[1] = -xi * (1 - zeta);           dN_zeta[1] = -xi * (1 - eta);
            dN_xi[2] = eta * (1 - zeta);            dN_eta[2] = xi * (1 - zeta);            dN_zeta[2] = -xi * eta;
            dN_xi[3] = -eta * (1 - zeta);           dN_eta[3] = (1 - xi) * (1 - zeta);      dN_zeta[3] = - (1 - xi) * eta;
            dN_xi[4] = - (1 - eta) * zeta;          dN_eta[4] = - (1 - xi) * zeta;          dN_zeta[4] = (1 - xi) * (1 - eta);
            dN_xi[5] = (1 - eta) * zeta;            dN_eta[5] = -xi * zeta;                 dN_zeta[5] = xi * (1 - eta);
            dN_xi[6] = eta * zeta;                  dN_eta[6] = xi * zeta;                  dN_zeta[6] = xi * eta;
            dN_xi[7] = -eta * zeta;                 dN_eta[7] = (1 - xi) * zeta;            dN_zeta[7] = (1 - xi) * eta;
            break;

        case GMSH_WEDGE:
            dN_xi[0] = -(1 - zeta);     dN_eta[0] = -(1 - zeta);    dN_zeta[0] = -(1 - xi - eta);
            dN_xi[1] = (1 - zeta);      dN_eta[1] = 0.0;            dN_zeta[1] = -xi;
            dN_xi[2] = 0.0;             dN_eta[2] = (1 - zeta);     dN_zeta[2] = -eta;
            dN_xi[3] = -zeta;           dN_eta[3] = -zeta;          dN_zeta[3] = (1 - xi - eta);
            dN_xi[4] = zeta;            dN_eta[4] = 0.0;            dN_zeta[4] = xi;
            dN_xi[5] = 0.0;             dN_eta[5] = zeta;           dN_zeta[5] = eta;
            break;

        case GMSH_PYRAMID:
            dN_xi[0] = 0.25 * ( - 1.0 + eta / (1 - zeta) );
            dN_xi[1] = 0.25 * ( 1.0 - eta / (1 - zeta) );
            dN_xi[2] = 0.25 * ( 1.0 + eta / (1 - zeta) );
            dN_xi[3] = 0.25 * ( - 1.0 - eta / (1 - zeta) );
            dN_xi[4] = 0.0;

            dN_eta[0] = 0.25 * ( - 1.0 + xi / (1 - zeta) );
            dN_eta[1] = 0.25 * ( - 1.0 - xi / (1 - zeta) );
            dN_eta[2] = 0.25 * ( 1.0 + xi / (1 - zeta) );
            dN_eta[3] = 0.25 * ( 1.0 - xi / (1 - zeta) );
            dN_eta[4] = 0.0;

            dN_zeta[0] = 0.25 * ( - 1.0 + (xi * eta) / pow(1 - zeta, 2) );
            dN_zeta[1] = 0.25 * ( - 1.0 - (xi * eta) / pow(1 - zeta, 2) );
            dN_zeta[2] = 0.25 * ( - 1.0 + (xi * eta) / pow(1 - zeta, 2) );
            dN_zeta[3] = 0.25 * ( - 1.0 - (xi * eta) / pow(1 - zeta, 2) );
            dN_zeta[4] = 1.0;
            break;

        default:
            fprintf(stderr, "Unknown element type: %d\n", elementType);
            exit(EXIT_FAILURE);
    }
}
/* ------------------------------------------------------------------------------------------------- */
int GetNumBasisFunctions(int elementType) {
    // printf("GetNumBasisFunctions called with elementType = %d\n", elementType);
    switch (elementType) {
        case GMSH_SEGMENT: return 2;      // 2 nodes
        case GMSH_TRIANGLE: return 3;     // 3 nodes
        case GMSH_QUADRANGLE: return 4;   // 4 nodes
        case GMSH_TETRAHEDRON: return 4;  // 4 nodes
        case GMSH_HEXAHEDRON: return 8;   // 8 nodes
        case GMSH_WEDGE: return 6;        // 6 nodes
        case GMSH_PYRAMID: return 5;      // 5 nodes
        default:
            fprintf(stderr, "Unknown element type %d (GetNumBasisFunctions)\n", elementType);
            exit(EXIT_FAILURE);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void ComputeGradientMatrixG(gelement3D *K, double G[3][8], double xi, double eta, double zeta) {

    double gradN_ref[3][8] = {0};
    GetReferenceGradients(K->type, gradN_ref[0], gradN_ref[1], gradN_ref[2], xi, eta, zeta);

    for (int j = 0; j < K->Nvertex; j++) { // Loop over nodes
        G[0][j] = K->JacobianInv[0][0] * gradN_ref[0][j] + K->JacobianInv[0][1] * gradN_ref[1][j] + K->JacobianInv[0][2] * gradN_ref[2][j];
        G[1][j] = K->JacobianInv[1][0] * gradN_ref[0][j] + K->JacobianInv[1][1] * gradN_ref[1][j] + K->JacobianInv[1][2] * gradN_ref[2][j];
        G[2][j] = K->JacobianInv[2][0] * gradN_ref[0][j] + K->JacobianInv[2][1] * gradN_ref[1][j] + K->JacobianInv[2][2] * gradN_ref[2][j];
    }

    // printf("Computed Gradients in Physical Space (G matrix):\n");
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < K->Nvertex; j++) {
    //         printf("G[%d][%d] = %f\n", i, j, G[i][j]);
    //     }
    // }
    // printf("\n");

    // printf("Computed invJ in Physical Space:\n");
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         printf("J_inv[%d][%d] = %f\n", i, j, K->JacobianInv[i][j]);
    //     }
    // }
    // printf("\n");

    // printf("Computed gradN_ref in Physical Space:\n");
    // for (int j = 0; j < K->Nvertex; j++) {
    //     for (int k = 0; k < 3; k++) {
    //         printf("gradN_ref[%d][%d] = %f\n", k, j, gradN_ref[k][j]);
    //     }
    // }
    // printf("\n");
}
/* ------------------------------------------------------------------------------------------------- */
int GetNumGaussPoints(int elementType) {
    switch (elementType) {
        case GMSH_TRIANGLE:      return 3;  // 1 Gauss point
        case GMSH_QUADRANGLE:    return 4;  // 2x2 Gauss quadrature
        case GMSH_TETRAHEDRON:   return 4;  // 1 centroid
        case GMSH_HEXAHEDRON:    return 8;  // 2x2x2 quadrature
        case GMSH_WEDGE:         return 6;  // 3 points (1 in base x 3 in height)
        case GMSH_PYRAMID:       return 5;  // 5 points
        default:
            fprintf(stderr, "Unknown element type %d (GetNumGaussPoint)\n", elementType);
            exit(EXIT_FAILURE);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void GetGaussPoint(int elementType, int gp, double *xi, double *eta, double *zeta, double *weight) {
    switch (elementType) {
        case GMSH_TRIANGLE: // 3-point quadrature - ok
            {
                double gp_points[3][2] = {
                    {0.5, 0.5},
                    {0.5, 0.0},
                    {0.0, 0.5}
                };
                double gp_weights[3] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};

                *xi = gp_points[gp][0];
                *eta = gp_points[gp][1];
                *zeta = 0.0;
                *weight = gp_weights[gp];
            }
            break;

        case GMSH_QUADRANGLE: // 2x2 quadrature -skipped, no modificaiton right now
            {
                double gp_points[4][2] = {
                    {-0.5773502692, -0.5773502692},
                    { 0.5773502692, -0.5773502692},
                    { 0.5773502692,  0.5773502692},
                    {-0.5773502692,  0.5773502692}
                };
                double gp_weights[4] = {1.0, 1.0, 1.0, 1.0};

                *xi = gp_points[gp][0];
                *eta = gp_points[gp][1];
                *zeta = 0.0;
                *weight = gp_weights[gp];
            }
            break;

        case GMSH_TETRAHEDRON: // 4-point quadrature - ok
            {
                double gp_points[4][3] = {
                    {0.13819660115, 0.13819660115, 0.13819660115},
                    {0.5854101966, 0.13819660115, 0.13819660115},
                    {0.13819660115, 0.5854101966, 0.13819660115},
                    {0.13819660115, 0.13819660115, 0.5854101966}
                };                
                double gp_weights[4] = {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0};

                *xi = gp_points[gp][0];
                *eta = gp_points[gp][1];
                *zeta = gp_points[gp][2];
                *weight = gp_weights[gp];
            }
            break;

        case GMSH_HEXAHEDRON: // 8-point quadrature (2x2x2)
            {
                double gp_points[8][3] = {
                    {-0.5773502692, -0.5773502692, -0.5773502692},
                    { 0.5773502692, -0.5773502692, -0.5773502692},
                    { 0.5773502692,  0.5773502692, -0.5773502692},
                    {-0.5773502692,  0.5773502692, -0.5773502692},
                    {-0.5773502692, -0.5773502692,  0.5773502692},
                    { 0.5773502692, -0.5773502692,  0.5773502692},
                    { 0.5773502692,  0.5773502692,  0.5773502692},
                    {-0.5773502692,  0.5773502692,  0.5773502692}
                };
                double gp_weights[8] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

                *xi = gp_points[gp][0];
                *eta = gp_points[gp][1];
                *zeta = gp_points[gp][2];
                *weight = gp_weights[gp];
            }
            break;

        case GMSH_WEDGE: // 6-point quadrature
            {
                double gp_points[6][3] = {
                    {0.5, 0.0, -0.5773502692},
                    {0.5, 0.5, -0.5773502692},
                    {0.0, 0.5, -0.5773502692},
                    {0.5, 0.0,  0.5773502692},
                    {0.5, 0.5,  0.5773502692},
                    {0.0, 0.5,  0.5773502692}
                };
                double gp_weights[6] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};

                *xi = gp_points[gp][0];
                *eta = gp_points[gp][1];
                *zeta = gp_points[gp][2];
                *weight = gp_weights[gp];
            }
            break;

        case GMSH_PYRAMID: // 5-point quadrature
            {
                double gp_points[5][3] = {
                    {0.0, 0.0, 0.25},
                    {0.5, 0.0, 0.25},
                    {0.0, 0.5, 0.25},
                    {0.5, 0.5, 0.25},
                    {0.25, 0.25, 0.75}
                };
                double gp_weights[5] = {0.1, 0.1, 0.1, 0.1, 0.6};

                *xi = gp_points[gp][0];
                *eta = gp_points[gp][1];
                *zeta = gp_points[gp][2];
                *weight = gp_weights[gp];
            }
            break;

        default:
            fprintf(stderr, "Unknown element type %d\n", elementType);
            exit(EXIT_FAILURE);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void ComputePhysicalCoordinates(double *x, double *y, double *z, gelement3D *K, double xi, double eta, double zeta) {
    int nBaseNodes = K->Nvertex; // Get the number of basis functions for the element type
    double N[8] = {0};           // Maximum of 8 nodes (e.g., hexahedron)

    ComputeBasisFunctions(K->type, N, xi, eta, zeta, nBaseNodes);

    *x = 0.0;
    *y = 0.0;
    *z = 0.0;

    // Accumulate physical coordinates using the shape functions
    for (int n = 0; n < nBaseNodes; n++) {
        *x += N[n] * K->Ver[n].x;
        *y += N[n] * K->Ver[n].y;
        *z += N[n] * K->Ver[n].z;
    }
}
/* ------------------------------------------------------------------------------------------------- */
//                                POISSON EQUATION CACLCULATION
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Poisson(mesh *this, double f(double, double), triplet *p, double *b, double *dirichletValues, double (*psi)(double, double, double, double), double (*ur)(double, double, double, double), IntegrationType integrationType)
{
    int i, j, k;
    memset(b, 0, this->NPoints * sizeof(double));

    for (k = 0; k < this->NElements; k++) {
        
        gelement K = Gelement(this, k);
        double G[2][3];
        computeGradient(&K, G);

        for (i = 0; i < 3; i++) {
            double integralValue = CalculateIntegral(this, K, f, integrationType, i);
            for (j = 0; j < 3; j++) {
                
                double qval = K.dV * (G[0][i] * G[0][j] + G[1][i] * G[1][j]);
                Triplet_Add(p, K.idx[i], K.idx[j], qval);
            }
            b[K.idx[i]] += integralValue;
        }
    }

#ifdef RHS_INTEGRAL
    printf("Load vector before boundary conditon:\n");
    for ( i = 0; i < this->NPoints; i++)
    {
        printf("b[%d] = %f\n",i,b[i]);
    }
#endif

//-------------------------------------------------------------------------------------------
    bool hasRobinBC = false;
    double *ContributionsROB = NULL;
    
    for (i = 0; i < this->NPoints; i++) {
        if (this->PointMark[i] < ROBIN_BC && this->PointMark[i] > NEUMANN_BC) {
            hasRobinBC = true;
            break;
        }
    }

    bool hasNeumannBC = false;
    double *ContributionsNEU = NULL;

    for (i = 0; i < this->NPoints; i++) {
        if (this->PointMark[i] < NEUMANN_BC && this->PointMark[i] > DIRICHLET_BC) {
            hasNeumannBC = true;
            break;
        }
    }

//-------------------------------------------------------------------------------------------

// Robin boundary condion

if (hasRobinBC) {
    ContributionsROB = calloc(this->NPoints, sizeof(double));
    for ( i = 0; i < this->NBEdges; i++) {
        Mesh_GetBndrSide S = MeshGetBoundarySide(this, i);
        if (this->PointMark[S.idxA] < ROBIN_BC || this->PointMark[S.idxB] < ROBIN_BC) {
            
            double urValueA = ur(this->x[S.idxA], this->y[S.idxA], S.nn[0], S.nn[1]);
            double urValueB = ur(this->x[S.idxB], this->y[S.idxB], S.nn[0], S.nn[1]);
            double matrixValue = S.ds * alfa / 2.0;

            // Load vector contributions
            ContributionsROB[S.idxA] += S.ds * alfa * urValueA / 2.0;
            ContributionsROB[S.idxB] += S.ds * alfa * urValueB / 2.0;

            #ifdef ROBIN_VALUES
            printf("\n");
            printf("idxA: %d\n",S.idxA);
            printf("idxB: %d\n",S.idxB);
            printf("idxA__con1: %f\n",S.ds * alfa * urValueA / 2.0);
            printf("idxB__con1: %f\n",S.ds * alfa * urValueB / 2.0);
            printf("\n");
            #endif
                
            // Update the stiffness matrix
            Triplet_Add(p, S.idxA, S.idxA, matrixValue);
            Triplet_Add(p, S.idxB, S.idxB, matrixValue);

        }
    }
} 

// Neumann boundary conditon

if (hasNeumannBC) {
    ContributionsNEU = calloc(this->NPoints, sizeof(double));
    for ( i = 0; i < this->NBEdges; i++) {
        Mesh_GetBndrSide S = MeshGetBoundarySide(this, i);
        if (this->PointMark[S.idxA] < NEUMANN_BC || this->PointMark[S.idxB] < NEUMANN_BC) {
            
            double psiValueA = psi(this->x[S.idxA], this->y[S.idxA], S.nn[0], S.nn[1]);
            double psiValueB = psi(this->x[S.idxB], this->y[S.idxB], S.nn[0], S.nn[1]);

            ContributionsNEU[S.idxA] += S.ds * psiValueA / 2.0;
            ContributionsNEU[S.idxB] += S.ds * psiValueB / 2.0;

        }
    }
}
    
for ( i = 0; i < this->NPoints; i++) {
    if (hasRobinBC) {    
        if (this->PointMark[i] > NEUMANN_BC ) {
            b[i] += ContributionsROB[i];
        }
    }
    if (hasNeumannBC) { 
        if (this->PointMark[i] > DIRICHLET_BC && this->PointMark[i] < NEUMANN_BC) {
            b[i] += ContributionsNEU[i];
        }
    }
}

if (hasRobinBC) {
    free(ContributionsROB);
}
if (hasNeumannBC) {
    free(ContributionsNEU);
}
//-------------------------------------------------------------------------------------------

//Dirichlet boundary conditions
    for (i = 0; i < this->NPoints; i++) {
        if (this->PointMark[i] == DIRICHLET_BC) {
            for (int nz = 0; nz < p->nz; nz++) {
                if (p->tripletI[nz] == i) {
                    int col = p->J[nz];
                    if (col != i) {
                        b[col] -= p->VAL[nz] * dirichletValues[i];
                    }
                }
            }
            
            // Set the diagonal entry to 1 and off-diagonal entries to 0
            for (j = 0; j < p->nz; j++) {
                if (p->tripletI[j] == i || p->J[j] == i) {
                    p->VAL[j] = 0.0;
                }
            }
            Triplet_Add(p, i, i, 1.0); // Diagonal is set to 1
            b[i] = dirichletValues[i]; 
        }
    }

    #ifdef LOAD_VECTOR
    printf("Load Vector b Final:\n");
    for (i = 0; i < this->NPoints; i++) {
        printf("b[%d] = %f\n", i, b[i]);
    }
    #endif

//-------------------------------------------------------------------------------------------

    Triplet_Sort(p, 0, p->nz - 1);
    Triplet_Unique(p);

//-------------------------------------------------------------------------------------------

    // Write the matrix to a file
    // #define MATRIX_CHECK
    #ifdef MATRIX_CHECK
    FILE *fid = fopen("TEST_POISSON.dat", "w+");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file to write matrix\n");
        return;
    }
    for (i = 0; i < p->nz; i++) {
        fprintf(fid, "%10d %10d %24.15g\n", p->I[i] + 1, p->J[i] + 1, p->VAL[i]);
    }
    fclose(fid);
    #endif

}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Poisson_3D(mesh *this, double (*f)(double, double, double), triplet *p, double *b, double *dirichletValues, double (*psi)(double, double, double, double, double, double), double (*ur)(double, double, double, double, double, double)) {
    
    int i, j, k;
    memset(b, 0, this->NPoints * sizeof(double));
    
    gelement3D K;
    
    double G[3][8];                 // Gradient matrix G for 3D elements
    double N[8];                    // Shape functions (max 8 nodes per element)
    double xi, eta, zeta, weight;   // Gauss points and weight
    
    for ( k = this->NTriangles + this->NQuads; k < this->NElements; k++ ) {
        
        GetElement3D(this, &K, k);
       
        int nBaseNodes = GetNumBasisFunctions(K.type);
        int nGaussPoints = GetNumGaussPoints(K.type);

        for (int i = 0; i < nBaseNodes; i++)
        {
            double val = 0.0;
            for (int gp = 0; gp < nGaussPoints; gp++)
            {
                GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);

                double x, y, z;
                ComputePhysicalCoordinates(&x, &y, &z, &K, xi, eta, zeta);
                ComputeBasisFunctions(K.type, N, xi, eta, zeta, nBaseNodes);
                computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, K.Jacobian, K.JacobianInv, &K.detJ, G, xi, eta, zeta);

                val += K.detJ * weight * f(x,y,z) * N[i];
            }
            b[K.idxNode[i]] += val;
        }

        for ( i = 0; i < nBaseNodes; i++)
        {
            for ( j = 0; j < nBaseNodes; j++)
            {
                double stiff_val = 0.0;
                for (int gp = 0; gp < nGaussPoints; gp++)
                {
                    GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);
                    computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, K.Jacobian, K.JacobianInv, &K.detJ, G, xi, eta, zeta);
                    double dotproduct = G[0][i] * G[0][j] + G[1][i] * G[1][j] + G[2][i] * G[2][j];
                    stiff_val += K.detJ * weight * dotproduct;
                }
                Triplet_Add(p, K.idxNode[i], K.idxNode[j], stiff_val);
            }
        }
    }

    //-------------------------------------------------------------------------------------------
    // bool hasRobinBC = false, hasNeumannBC = false;
    // double *ContributionsROB = NULL, *ContributionsNEU = NULL;

    // for (int i = 0; i < this->NTriangles + this->NQuads; i++) {
    //     if (this->ElementMark[i] < ROBIN_BC && this->ElementMark[i] > NEUMANN_BC) {
    //         hasRobinBC = true;
    //         break;
    //     }
    // }

    // for (int i = 0; i < this->NTriangles + this->NQuads; i++) {
    //     if (this->ElementMark[i] < NEUMANN_BC && this->ElementMark[i] > DIRICHLET_BC) {
    //         hasNeumannBC = true;
    //         break;
    //     }
    // }
    //-------------------------------------------------------------------------------------------
    #define NEU
    #ifdef NEU
    // Neumann boundary conditions
    // if (hasNeumannBC) {
        for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
            if (this->ElementMark[elem] < NEUMANN_BC && this->ElementMark[elem] > DIRICHLET_BC) {
                GetElement3D(this, &K, elem);
                // printf("Element index: %d\n", elem);
                if (K.type != GMSH_TRIANGLE && K.type != GMSH_QUADRANGLE) {
                    fprintf(stderr, "Unsupported face type %d in Neumann BC\n", K.type);
                    continue;
                }

                int nVertices = GetNumBasisFunctions(K.type);
                int nGaussPoints = GetNumGaussPoints(K.type);

                for (int i = 0; i < nVertices; i++) {
                    double uN = 0.0;
                    int node = K.idxNode[i];

                    for (int gp = 0; gp < nGaussPoints; gp++) {
                        GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);

                        double gradN_xi[8] = {0.0}, gradN_eta[8] = {0.0};
                        GetReferenceGradients(K.type, gradN_xi, gradN_eta, NULL, xi, eta, zeta);

                        double J_face = 0.0, normal_gp[3] = {0.0};
                        computeSurfaceJacobianAndNormal(K.Ver, nVertices, gradN_xi, gradN_eta, &J_face, normal_gp, NULL);

                        double x, y, z;
                        ComputePhysicalCoordinates(&x, &y, &z, &K, xi, eta, zeta);
                        ComputeBasisFunctions(K.type, N, xi, eta, zeta, nVertices);

                        double psiValue = psi(x, y, z, normal_gp[0], normal_gp[1], normal_gp[2]);
                        double incr = J_face * weight * N[i] * psiValue;
                        uN += incr;

                        // printf("gp = %d | i = %d | J_face = %g | w = %g | N = %g | psi = %g | incr = %g\n", gp, i, J_face, weight, N[i], psiValue, incr);
                    }

                    b[node] += uN;
                    // printf("node = %d | b[node] =  %g\n", node, b[node]);
                }
                // printf("\n");
            }
        }
    // }
    #endif

    // printf("Neumann:\n");
    // for (int i = 0; i < this->NPoints; i++)
    // {
    //     printf("b[%d] = %g\n", i, b[i]);
    // }
    // printf("\n");
    // #define ROB
    #ifdef ROB
    // Robin boundary conditions
    if (hasRobinBC) {
        ContributionsROB = calloc(this->NPoints, sizeof(double));
    
        for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
            if (this->ElementMark[elem] < ROBIN_BC && this->ElementMark[elem] > NEUMANN_BC) {
    
                GetElement3D(this, &K, elem);
    
                if (K.type != GMSH_TRIANGLE && K.type != GMSH_QUADRANGLE) {
                    fprintf(stderr, "Invalid element type %d for Robin BC\n", K.type);
                    continue;
                }
    
                int nVertices = GetNumBasisFunctions(K.type);
                int nGaussPoints = GetNumGaussPoints(K.type);
    
                for (int gp = 0; gp < nGaussPoints; gp++) {
                    GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);
                    computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, K.Jacobian, K.JacobianInv, &K.detJ, G, xi, eta, zeta);
    
                    double x, y, z;
                    ComputePhysicalCoordinates(&x, &y, &z, &K, xi, eta, zeta);
                    ComputeBasisFunctions(K.type, N, xi, eta, zeta, nVertices);
    
                    double urValue = ur(x, y, z, K.normal[0], K.normal[1], K.normal[2]);
    
                    for (int i = 0; i < nVertices; i++) {
                        int node_i = K.idxNode[i];
    
                        // Load vector contribution
                        double load_val = K.area * weight * alfa * urValue * N[i];
                        // ContributionsROB[node_i] += load_val;
    
                        for (int j = 0; j < nVertices; j++) {
                            int node_j = K.idxNode[j];
    
                            // Stiffness matrix contribution
                            double stiff_val = K.area * weight * alfa * N[i] * N[j];
                            Triplet_Add(p, node_i, node_j, stiff_val);
                        }
                    }
                }
            }
        }
    }
    #endif

    // if (hasRobinBC) {
        // ContributionsROB = calloc(this->NPoints, sizeof(double));
        #define ROBIN_BETTER
        #ifdef ROBIN_BETTER
        for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
            if (this->ElementMark[elem] < ROBIN_BC && this->ElementMark[elem] > NEUMANN_BC) {
    
                GetElement3D(this, &K, elem);
    
                if (K.type != GMSH_TRIANGLE && K.type != GMSH_QUADRANGLE) {
                    fprintf(stderr, "Invalid element type %d for Robin BC\n", K.type);
                    continue;
                }
    
                int nVertices = GetNumBasisFunctions(K.type);
                int nGaussPoints = GetNumGaussPoints(K.type);
    
                for (int i = 0; i < nVertices; i++) {
                    int node_i = K.idxNode[i];
                    double robinLoad = 0.0;
    
                    for (int gp = 0; gp < nGaussPoints; gp++) {
                        GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);
    
                        double gradN_xi[8] = {0.0}, gradN_eta[8] = {0.0};
                        GetReferenceGradients(K.type, gradN_xi, gradN_eta, NULL, xi, eta, zeta);
    
                        double J_face = 0.0, normal_gp[3] = {0.0};
                        computeSurfaceJacobianAndNormal(K.Ver, nVertices, gradN_xi, gradN_eta, &J_face, normal_gp, NULL);
    
                        double x, y, z;
                        ComputePhysicalCoordinates(&x, &y, &z, &K, xi, eta, zeta);
                        ComputeBasisFunctions(K.type, N, xi, eta, zeta, nVertices);
    
                        // Robin condition with ur(x,y,z) at Gauss point
                        double urValue = ur(x, y, z, normal_gp[0], normal_gp[1], normal_gp[2]);
    
                        // Add Robin load contribution like Neumann
                        robinLoad += J_face * weight * alfa * urValue * N[i];
                    }
    
                    b[node_i] += robinLoad;
                }
    
                // Robin stiffness matrix contribution (N_i * N_j terms)
                for (int i = 0; i < nVertices; i++) {
                    int node_i = K.idxNode[i];
    
                    for (int j = 0; j < nVertices; j++) {
                        int node_j = K.idxNode[j];
                        double robinStiff = 0.0;
    
                        for (int gp = 0; gp < nGaussPoints; gp++) {
                            GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);
    
                            double gradN_xi[8] = {0.0}, gradN_eta[8] = {0.0};
                            GetReferenceGradients(K.type, gradN_xi, gradN_eta, NULL, xi, eta, zeta);
    
                            double J_face = 0.0, normal_gp[3] = {0.0};
                            computeSurfaceJacobianAndNormal(K.Ver, nVertices, gradN_xi, gradN_eta, &J_face, normal_gp, NULL);
    
                            ComputeBasisFunctions(K.type, N, xi, eta, zeta, nVertices);
    
                            // Robin weak contribution
                            robinStiff += J_face * weight * alfa * N[i] * N[j];
                        }
    
                        Triplet_Add(p, node_i, node_j, robinStiff);
                    }
                }
            }
        }
    // }
        #endif
    // Apply contributions from Robin and Neumann conditions to the load vector
    // for (int i = 0; i < this->NPoints; i++) {
    //     if (hasRobinBC && this->PointMark[i] > NEUMANN_BC) {
    //         b[i] += ContributionsROB[i];
    //     }
    //     if (hasNeumannBC && this->PointMark[i] > DIRICHLET_BC && this->PointMark[i] < NEUMANN_BC) {
    //         b[i] += ContributionsNEU[i];
    //     }
    // }
    // printf("Robin:\n");
    // for (int i = 0; i < this->NPoints; i++)
    // {
    //     printf("b[%d] = %g\n", i, b[i]);
    // }
    // printf("\n");

    // Free memory for boundary contributions
    // if (hasRobinBC) {
    //     free(ContributionsROB);
    // }
    // if (hasNeumannBC) {
    //     free(ContributionsNEU);
    // }
    // -------------------------------------------------------------------------------------------
    // Dirichlet boundary conditions
    // printf("Nelements2D: %d\n", this->NTriangles + this->NQuads);
    #define DIR
    #ifdef DIR
    for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
        if (this->ElementMark[elem] == DIRICHLET_BC) {

            GetElement3D(this, &K, elem);
            int nVertices = GetNumBasisFunctions(K.type);

            for (int i = 0; i < nVertices; i++)
            {
                int node = this->ElementListOfVertices[nVertices * elem + i];
                    for (int nz = 0; nz < p->nz; nz++)
                    {
                        if (p->tripletI[nz] == node) {
                            int col = p->J[nz];
                            if (col != node) {
                                b[col] -= p->VAL[nz] * dirichletValues[node];
                            }   
                        } 
                    }

                    // Set the diagonal entry to 1 and off-diagonal entries to 0
                    for (int j = 0; j < p->nz; j++) {
                        if (p->tripletI[j] == node || p->J[j] == node) {
                            p->VAL[j] = 0.0;
                        }
                    }
                    Triplet_Add(p, node, node, 1.0); // Diagonal is set to 1
                    b[node] = dirichletValues[node]; // Set the Dirichlet value in the load vector
            }
        }
    }
    
    #endif
    // #define LOAD_VECTOR3D
    #ifdef LOAD_VECTOR3D
    printf("Load Vector b Final:\n");
    for (i = 0; i < this->NPoints; i++) {
        printf("b[%d] = %f\n", i, b[i]);
    }
    // saveVectorToFile(b, this->NPoints, 1, "Load_Vector3D.txt");
    #endif
    //-------------------------------------------------------------------------------------------
    Triplet_Sort(p, 0, p->nz - 1);
    Triplet_Unique(p);

    #define MATRIX3D_CHECK
    #ifdef MATRIX3D_CHECK
    FILE *fid = fopen("TEST3D_POISSON.dat", "w+");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file to write matrix\n");
        return;
    }
    for (i = 0; i < p->nz; i++) {
        fprintf(fid, "%10d %10d %24.15g\n", p->tripletI[i] + 1, p->J[i] + 1, p->VAL[i]);
    }
    fclose(fid);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
//                                      NORMS CALCULATION
/* ------------------------------------------------------------------------------------------------- */
void computeGradientAtCentroid(mesh *M, double *solution, int elementIndex, double *gradX, double *gradY)
{
    gelement K = Gelement(M, elementIndex);

    double uA = solution[K.idx[0]];
    double uB = solution[K.idx[1]];
    double uC = solution[K.idx[2]];

    *gradX = ((K.B[1] - K.C[1]) * uA + (K.C[1] - K.A[1]) * uB + (K.A[1] - K.B[1]) * uC) / (2.0 * K.dV);
    *gradY = ((K.C[0] - K.B[0]) * uA + (K.A[0] - K.C[0]) * uB + (K.B[0] - K.A[0]) * uC) / (2.0 * K.dV);
}
/* ------------------------------------------------------------------------------------------------- */
double computeL2Error(mesh *M, double *numericalSolution, double (*analyticalSolution)(double, double))
{
    double errorSquaredSum = 0.0;

    for (int i = 0; i < M->NElements; i++) {
        gelement K = Gelement(M, i);

        double x_centroid = (K.A[0] + K.B[0] + K.C[0]) / 3.0;
        double y_centroid = (K.A[1] + K.B[1] + K.C[1]) / 3.0;

        // printf("(X,Y) = ( %g, %g )\n",x_centroid, y_centroid);

        double u_centroid = (numericalSolution[K.idx[0]] + numericalSolution[K.idx[1]] + numericalSolution[K.idx[2]]) / 3.0;

        double analytical_centroid = analyticalSolution(x_centroid, y_centroid);

        double errorSquared = pow(u_centroid - analytical_centroid, 2);

        errorSquaredSum += errorSquared * K.dV;
    }

    return sqrt(errorSquaredSum);
}
/* ------------------------------------------------------------------------------------------------- */
double computeH1SeminormError(mesh *M, double *numericalSolution, double (*AnalyticalSolution)(double, double), void (*CalculateAnalyticalGradient)(double, double, double*, double*))
{
    double l2Error = computeL2Error(M, numericalSolution, AnalyticalSolution);
    double gradErrorSquaredSum = 0.0;

    for (int i = 0; i < M->NElements; i++) {
        double gradX, gradY;
        computeGradientAtCentroid(M, numericalSolution, i, &gradX, &gradY);

        // Compute centroid
        gelement K = Gelement(M, i);
        double x_centroid = (K.A[0] + K.B[0] + K.C[0]) / 3.0;
        double y_centroid = (K.A[1] + K.B[1] + K.C[1]) / 3.0;

        // Analytical gradient
        double analGradX, analGradY;
        CalculateAnalyticalGradient(x_centroid, y_centroid, &analGradX, &analGradY);

        // Gradient error
        double gradErrorX = gradX - analGradX;
        double gradErrorY = gradY - analGradY;
        gradErrorSquaredSum += (gradErrorX * gradErrorX + gradErrorY * gradErrorY) * K.dV;
    }

    return sqrt(l2Error * l2Error + gradErrorSquaredSum);
}
/* ------------------------------------------------------------------------------------------------- */
void computeGradientAtCentroid3D(mesh *M, double *solution, int elementIndex, double *gradX, double *gradY, double *gradZ)
{
    gelement3D K;
    GetElement3D(M, &K, elementIndex); // Populate gelement3D structure for the given element index.

    // Extract solution values at the vertices
    double uA = solution[K.idxNode[0]];
    double uB = solution[K.idxNode[1]];
    double uC = solution[K.idxNode[2]];
    double uD = solution[K.idxNode[3]];

    // Compute the Jacobian determinant
    double J[3][3] = {
        {K.Ver[1].x - K.Ver[0].x, K.Ver[2].x - K.Ver[0].x, K.Ver[3].x - K.Ver[0].x},
        {K.Ver[1].y - K.Ver[0].y, K.Ver[2].y - K.Ver[0].y, K.Ver[3].y - K.Ver[0].y},
        {K.Ver[1].z - K.Ver[0].z, K.Ver[2].z - K.Ver[0].z, K.Ver[3].z - K.Ver[0].z},
    };

    double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    if (fabs(detJ) < 1e-12) {
        fprintf(stderr, "Error: Degenerate element with zero or near-zero determinant.\n");
        exit(EXIT_FAILURE);
    }

    double invDetJ = 1.0 / detJ;

    // Compute gradients of the shape functions
    double gradN[4][3]; // Gradients of shape functions at centroid
    gradN[0][0] = invDetJ * ((K.Ver[2].y - K.Ver[3].y) * (K.Ver[1].z - K.Ver[3].z) - (K.Ver[2].z - K.Ver[3].z) * (K.Ver[1].y - K.Ver[3].y));
    gradN[0][1] = invDetJ * ((K.Ver[2].z - K.Ver[3].z) * (K.Ver[1].x - K.Ver[3].x) - (K.Ver[2].x - K.Ver[3].x) * (K.Ver[1].z - K.Ver[3].z));
    gradN[0][2] = invDetJ * ((K.Ver[2].x - K.Ver[3].x) * (K.Ver[1].y - K.Ver[3].y) - (K.Ver[2].y - K.Ver[3].y) * (K.Ver[1].x - K.Ver[3].x));

    gradN[1][0] = invDetJ * ((K.Ver[3].y - K.Ver[2].y) * (K.Ver[0].z - K.Ver[2].z) - (K.Ver[3].z - K.Ver[2].z) * (K.Ver[0].y - K.Ver[2].y));
    gradN[1][1] = invDetJ * ((K.Ver[3].z - K.Ver[2].z) * (K.Ver[0].x - K.Ver[2].x) - (K.Ver[3].x - K.Ver[2].x) * (K.Ver[0].z - K.Ver[2].z));
    gradN[1][2] = invDetJ * ((K.Ver[3].x - K.Ver[2].x) * (K.Ver[0].y - K.Ver[2].y) - (K.Ver[3].y - K.Ver[2].y) * (K.Ver[0].x - K.Ver[2].x));

    gradN[2][0] = invDetJ * ((K.Ver[0].y - K.Ver[3].y) * (K.Ver[1].z - K.Ver[3].z) - (K.Ver[0].z - K.Ver[3].z) * (K.Ver[1].y - K.Ver[3].y));
    gradN[2][1] = invDetJ * ((K.Ver[0].z - K.Ver[3].z) * (K.Ver[1].x - K.Ver[3].x) - (K.Ver[0].x - K.Ver[3].x) * (K.Ver[1].z - K.Ver[3].z));
    gradN[2][2] = invDetJ * ((K.Ver[0].x - K.Ver[3].x) * (K.Ver[1].y - K.Ver[3].y) - (K.Ver[0].y - K.Ver[3].y) * (K.Ver[1].x - K.Ver[3].x));

    gradN[3][0] = invDetJ * ((K.Ver[1].y - K.Ver[2].y) * (K.Ver[0].z - K.Ver[2].z) - (K.Ver[1].z - K.Ver[2].z) * (K.Ver[0].y - K.Ver[2].y));
    gradN[3][1] = invDetJ * ((K.Ver[1].z - K.Ver[2].z) * (K.Ver[0].x - K.Ver[2].x) - (K.Ver[1].x - K.Ver[2].x) * (K.Ver[0].z - K.Ver[2].z));
    gradN[3][2] = invDetJ * ((K.Ver[1].x - K.Ver[2].x) * (K.Ver[0].y - K.Ver[2].y) - (K.Ver[1].y - K.Ver[2].y) * (K.Ver[0].x - K.Ver[2].x));

    // Compute the gradient of the solution
    *gradX = gradN[0][0] * uA + gradN[1][0] * uB + gradN[2][0] * uC + gradN[3][0] * uD;
    *gradY = gradN[0][1] * uA + gradN[1][1] * uB + gradN[2][1] * uC + gradN[3][1] * uD;
    *gradZ = gradN[0][2] * uA + gradN[1][2] * uB + gradN[2][2] * uC + gradN[3][2] * uD;
}
/* ------------------------------------------------------------------------------------------------- */
double computeL2Error3D(mesh *M, double *numericalSolution, double (*analyticalSolution)(double, double, double))
{
    double errorSquaredSum = 0.0;

    for (int i = M->NTriangles + M->NQuads; i < M->NElements; i++) {
        gelement3D K;
        GetElement3D(M, &K, i);

        // Calculate the centroid of the tetrahedral element
        double x_centroid = 0.0, y_centroid = 0.0, z_centroid = 0.0;
        
        for (int j = 0; j < K.Nvertex; j++) {
            x_centroid += K.Ver[j].x;
            y_centroid += K.Ver[j].y;
            z_centroid += K.Ver[j].z;
        }

        x_centroid /= K.Nvertex;
        y_centroid /= K.Nvertex;
        z_centroid /= K.Nvertex;

        // Interpolating the numerical solution at the centroid
        double u_centroid = 0.0;
        for (int j = 0; j < K.Nvertex; j++) {
            u_centroid += numericalSolution[K.idxNode[j]];
        }
        u_centroid /= K.Nvertex;

        // Analytical solution at the centroid
        double analytical_centroid = analyticalSolution(x_centroid, y_centroid, z_centroid);

        // Compute the squared error and accumulate
        double errorSquared = pow(u_centroid - analytical_centroid, 2);
        errorSquaredSum += errorSquared * K.volume;
    }

    return sqrt(errorSquaredSum);
}
/* ------------------------------------------------------------------------------------------------- */
double computeH1SeminormError3D(mesh *M, double *numericalSolution, double (*AnalyticalSolution)(double, double, double), void (*CalculateAnalyticalGradient)(double, double, double, double*, double*, double*))
{
    double l2Error = computeL2Error3D(M, numericalSolution, AnalyticalSolution);
    double gradErrorSquaredSum = 0.0;

    for (int i = 0; i < M->NElements; i++) {
        double gradX, gradY, gradZ;
        computeGradientAtCentroid3D(M, numericalSolution, i, &gradX, &gradY, &gradZ);

        gelement3D K;
        GetElement3D(M, &K, i);

        double x_centroid = (K.Ver[0].x + K.Ver[1].x + K.Ver[2].x + K.Ver[3].x) / 4.0;
        double y_centroid = (K.Ver[0].y + K.Ver[1].y + K.Ver[2].y + K.Ver[3].y) / 4.0;
        double z_centroid = (K.Ver[0].z + K.Ver[1].z + K.Ver[2].z + K.Ver[3].z) / 4.0;

        double analGradX, analGradY, analGradZ;
        CalculateAnalyticalGradient(x_centroid, y_centroid, z_centroid, &analGradX, &analGradY, &analGradZ);

        double gradErrorX = gradX - analGradX;
        double gradErrorY = gradY - analGradY;
        double gradErrorZ = gradZ - analGradZ;
        gradErrorSquaredSum += (gradErrorX * gradErrorX + gradErrorY * gradErrorY + gradErrorZ * gradErrorZ) * K.volume;
    }

    return sqrt(l2Error * l2Error + gradErrorSquaredSum);
}
/* ------------------------------------------------------------------------------------------------- */
//                                     OHYB NOSNIKU - 2D case
/* ------------------------------------------------------------------------------------------------- */
void initializeMaterialProperties(mesh *M, MaterialProperties *materials)
{
    if (!materials) return;

    int startIdx = (M->dim == 3) ? (M->NTriangles + M->NQuads) : 0;

    for (int i = startIdx; i < M->NElements; i++)
    {
        int marker = M->ElementMark[i];
        switch (marker) {
            case STEEL:
                materials[i].E = 210e9;
                materials[i].nu = 0.3;
                materials[i].rho = 7800;
                materials[i].model = SVK; // SVK

                // materials[i].E = 1.6e6;
                // materials[i].nu = 0.4;
                // materials[i].rho = 1250;
                // materials[i].model = SVK;// SVK My choosing

                // materials[i].E = 70e9;
                // materials[i].nu = 0.32;
                // materials[i].rho = 2710;
                // materials[i].model = SVK; //NEOHOOKE
                break;
            case ALUMINUM:
                materials[i].E = 70e9;
                materials[i].nu = 0.32;
                materials[i].rho = 2710;
                materials[i].model = SVK; //NEOHOOKE
                break;
            case PVC:
                materials[i].E = 3.2e9;
                materials[i].nu = 0.38;
                materials[i].rho = 1450;
                materials[i].model = NEOHOOKE; //NEOHOOKE
                break;
            case SILICONE:
                materials[i].E = 0.0255e9; //AZO MATERIALS SILICONE RUBBER
                materials[i].nu = 0.48;
                materials[i].rho = 1700;
                materials[i].model =  NEOHOOKE; //NEOHOOKE
                break;
            case MAHAGONY:
                // https://www.fpl.fs.usda.gov/documnts/fplgtr/fplgtr190/chapter_05.pdf

                materials[i].E1 = 7.9e9;
                materials[i].E2 = 0.8769e9;
                materials[i].E3 = 0.395e9;

                materials[i].nu12 = 0.297;
                materials[i].nu13 = 0.641;
                materials[i].nu23 = 0.604;
                
                materials[i].G12 = 0.6952e9;
                materials[i].G13 = 0.4661e9;
                materials[i].G23 = 0.1659e9;
                
                materials[i].rho = 650;
                materials[i].model = ORTHOTROPIC;
                break;
                default:
                printf("Warning: Unrecognized material marker %d for element %d. Using default material (STEEL) properties.\n", marker, i);
                exit(EXIT_FAILURE);
                break;
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void calculateBMatrix(const gelement *K, double B[3][6])
{
    
    double G[2][3];
    computeGradient(K, G);

    memset(B, 0, sizeof(double) * 3 * 6);
    
    // Populate the B matrix using the gradients of shape functions
    for (int i = 0; i < 3; i++) {
        
        B[0][2*i] = G[0][i];      // dN_i/dx
        B[1][2*i] = 0.0;
        B[2][2*i] = G[1][i];      // dN_i/dy for shear deformation
        
        B[0][2*i+1] = 0.0;
        B[1][2*i+1] = G[1][i];    // dN_i/dy
        B[2][2*i+1] = G[0][i];    // dN_i/dx for shear deformation
    }

    #ifdef Bmatrix
    printf("Element %d B Matrix:\n", K->idx[0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            printf("%12.5e ", B[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void computeElasticityMatrix(MaterialProperties mat, double D[3][3])
{
    double E = mat.E;
    double nu = mat.nu;
    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));

    D[0][0] = lambda + 2 * mu;
    D[0][1] = lambda;
    D[0][2] = 0.0;
    D[1][0] = lambda;
    D[1][1] = lambda + 2 * mu;
    D[1][2] = 0.0;
    D[2][0] = 0.0;
    D[2][1] = 0.0;
    D[2][2] = mu; // mu or 2*mu

    #ifdef Dmatrix
    printf("D Matrix:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%12.5e ", D[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void PrintTriplet(const char *label, triplet *this) {
    printf("\n%s:\n", label);
    printf("nz (non-zero count): %d\n", this->nz);
    for (int i = 0; i < this->nz; i++) {
        printf("I: %d, J: %d, VAL: %.10e\n", this->tripletI[i], this->J[i], this->VAL[i]);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void computeLocalStiffnessMatrix(mesh *this, int elementIndex, double localK[6][6], MaterialProperties *materials) {
    gelement K = Gelement(this, elementIndex);
    double B[3][6], D[3][3], G[2][3];
    MaterialProperties mat = materials[elementIndex];

    computeGradient(&K, G);
    calculateBMatrix(&K, B);
    computeElasticityMatrix(mat, D);

    memset(localK, 0, sizeof(double) * 36);  // Initialize local stiffness matrix
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    localK[i][j] += (B[a][i] * D[a][b] * B[b][j]) * K.dV;
                }
            }
        }
    }
    
    #ifdef LOC_MATRIX
        printf("Local stiffness matrix ke for element %d:\n", elementIndex);
        printf("Node Coordinates:\n");
        for (int idx = 0; idx < 3; idx++) {  // Assuming each element K has 3 nodes, typical for triangular elements
            printf("Node %d: (x = %f, y = %f)\n", K.idx[idx], this->x[K.idx[idx]], this->y[K.idx[idx]]);
        }
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                printf("%12.5e ", localK[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    #endif
    #ifdef SAVE_LOC
    // Open file in append mode to save each local matrix consecutively
    FILE *fid = fopen("Loc_matrix_all.dat", "a");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file to write matrix\n");
        return;
    }

    fprintf(fid, "Local stiffness matrix for element %d:\n", elementIndex);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            fprintf(fid, "%12.5e ", localK[i][j]);
        }
        fprintf(fid, "\n");
    }
    fprintf(fid, "\n");  // Add a newline between matrices
    fclose(fid);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void assembleGlobalStiffnessMatrix(mesh *this, triplet *p, MaterialProperties *materials) {
    double localK[6][6];

    for (int k = 0; k < this->NElements; k++) {
        gelement K = Gelement(this, k);
        computeLocalStiffnessMatrix(this, k, localK, materials);

        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                int global_i = 2 * K.idx[i / 2] + (i % 2); //(1 - i % 2)
                int global_j = 2 * K.idx[j / 2] + (j % 2); //(1 - j % 2)
                Triplet_Add(p, global_i, global_j, localK[i][j]); 
            }
        }
    }
    #ifdef MATRIX_NO_BC
    FILE *fid = fopen("Matrix_NO_BC.dat", "w+");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file to write matrix\n");
        return;
    }
    for (int i = 0; i < p->nz; i++) {
        fprintf(fid, "%10d %10d %24.15g\n", p->tripletI[i]+1, p->J[i]+1, p->VAL[i]);
    }
    fclose(fid);
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void applyGravity(mesh *this, double *b, MaterialProperties *materials) {
    double gravity = -9.81;
    for (int k = 0; k < this->NElements; k++) {
        gelement K = Gelement(this, k);
        MaterialProperties mat = materials[k];
        double elementWeight = mat.rho * K.dV * gravity;
        double loadPerNode = elementWeight / 3;

        for (int i = 0; i < 3; i++) {
            b[2 * K.idx[i] + 1] += loadPerNode;
        }
    }
    #ifdef LOAD_GRAVITY
    printf("After Gravitational load:\n");
    for (int i = 0; i < this->NPoints; i++)
    {
        // printf("b[%d] = %fb[%d] = %f\n",2*i,b[2*i],2*i+1,b[2*i+1]);
        printf("bx[%d] = %f\n",2*i,b[2*i]);
        printf("by[%d] = %f\n",2*i+1,b[2*i+1]);
    }
    #endif
}
/* ------------------------------------------------------------------------------------------------- */
void applyBoundaryConditions(mesh *this, double *b, double *dirichletValues, triplet *p, double (*psi_force)(double, double, double, double, double, double *, double *), double angle_degrees) 
{

    for (int i = 0; i < this->NBEdges; i++) {
        Mesh_GetBndrSide S = MeshGetBoundarySide(this, i);
        if (this->PointMark[S.idxA] == NEUMANN_LOAD || this->PointMark[S.idxB] == NEUMANN_LOAD) {
            double angle_rad = angle_degrees * (M_PI / 180.0);
            double cos_theta = cos(angle_rad);
            double sin_theta = sin(angle_rad);

            double force_x, force_y;
            double magnitude = psi_force(this->x[S.idxA], this->y[S.idxA], S.nn[0], S.nn[1], angle_degrees, &force_x, &force_y);
            
            double x_force = magnitude * (cos_theta * S.nn[1] - sin_theta * S.nn[0]) * S.ds;
            double y_force = magnitude * (cos_theta * S.nn[0] + sin_theta * S.nn[1]) * S.ds;

            b[2 * S.idxA] += x_force / 2;
            b[2 * S.idxA + 1] += y_force / 2;
            b[2 * S.idxB] += x_force / 2;
            b[2 * S.idxB + 1] += y_force / 2;

        }
        else if (this->PointMark[S.idxA] == NEUMANN_FREE || this->PointMark[S.idxB] == NEUMANN_FREE)
        {
            double x_force = 0.0;
            double y_force = 0.0;
            
            b[2 * S.idxA] += x_force / 2;
            b[2 * S.idxA + 1] += y_force / 2;
            b[2 * S.idxB] += x_force / 2;
            b[2 * S.idxB + 1] += y_force / 2;
        }
    }

    // Apply Dirichlet boundary conditions
    for (int i = 0; i < this->NPoints; i++) {
        if (this->PointMark[i] == DIRICHLET_BC) {
            int ui = 2 * i;
            int vi = 2 * i + 1;

            for (int nz = 0; nz < p->nz; nz++) {
                if (p->tripletI[nz] == ui || p->tripletI[nz] == vi || p->J[nz] == ui || p->J[nz] == vi) {
                    p->VAL[nz] = 0.0;
                }
            }

            // Add diagonal entries to enforce Dirichlet conditions
            Triplet_Add(p, ui, ui, 1.0);
            Triplet_Add(p, vi, vi, 1.0);

            // Set known Dirichlet displacements in the load vector
            b[ui] = dirichletValues[ui];
            b[vi] = dirichletValues[vi];
        }
    }

    Triplet_Sort(p, 0, p->nz - 1);
    Triplet_Unique(p);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Statics(mesh *this, double *b, triplet *p, double *dirichletValues, double (*psi_force)(double, double, double, double, double, double *, double *), double angle_degrees, int gravity) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    memset(b, 0, this->NPoints * sizeof(double) * 2);

    assembleGlobalStiffnessMatrix(this, p, materials);

    if (gravity)
    {
        applyGravity(this, b, materials);
    }
    
    applyBoundaryConditions(this, b, dirichletValues, p, psi_force, angle_degrees);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
//                                     DYNAMICS - 2D case
/* ------------------------------------------------------------------------------------------------- */
void computeConsistentMassMatrix(mesh *this, triplet *massTriplet) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    for (int e = 0; e < this->NElements; e++) {
        gelement K = Gelement(this, e);
        double localM[6][6] = {0};
        double rho = materials[e].rho;
        double multiplier = ( rho * K.dV ) / 12.0;

        // Compute local mass matrix entries
        /*
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                localM[i][j] = ((i % 2 == j % 2) && (i / 2 == j / 2) ? 2 : 1) * multiplier;
            }
        }
        */

        // Fill in the consistent mass matrix
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                localM[2 * i][2 * j] = (i == j ? 2 : 1) * multiplier;  // x-x components
                localM[2 * i + 1][2 * j + 1] = (i == j ? 2 : 1) * multiplier;  // y-y components
            }
        }

        #ifdef LOC_MASS
        printf("Local mass matrix for element %d:\n", e);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                printf("%12.5e ", localM[i][j]);
                // printf("%12.5e ", ( 1 / multiplier ) * localM[i][j]);s
            }
            printf("\n");
        }
        printf("\n");
        #endif        

        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                int global_i = 2 * K.idx[i / 2] + (i % 2);
                int global_j = 2 * K.idx[j / 2] + (j % 2);
                Triplet_Add(massTriplet, global_i, global_j, localM[i][j]);
            }
        }

        #ifdef MASS_MATRIX
            printf("Local mass matrix for element %d:\n", e);
            printf("Node Coordinates:\n");
            for (int idx = 0; idx < 3; idx++) {
                printf("Node %d: (x = %f, y = %f)\n", K.idx[idx], this->x[K.idx[idx]], this->y[K.idx[idx]]);
            }
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    printf("%12.5e ", localM[i][j]); //(12 / (K.dV * rho)) *
                }
                printf("\n");
            }
            printf("\n");
        #endif
        
    }

    Triplet_Sort(massTriplet, 0, massTriplet->nz - 1);
    Triplet_Unique(massTriplet);

    free(materials);

}
/* ------------------------------------------------------------------------------------------------- */
void assembleGlobalLumpedMassMatrix(mesh *this, triplet *massTriplet) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    for (int e = 0; e < this->NElements; e++) {
        gelement K = Gelement(this, e);
        double localM[6][6] = {0};
        double rho = materials[e].rho;
        double multiplier = (rho * K.dV) / 3.0;

        // Step 3: Fill the lumped mass matrix (diagonal only) - nodes of element only
        for (int i = 0; i < 3; i++) {

            localM[2 * i][2 * i] = multiplier;      // x-x component (node i)
            localM[2 * i + 1][2 * i + 1] = multiplier;  // y-y component (node i)
        }

        for (int i = 0; i < 6; i++) {
            int global_i = 2 * K.idx[i / 2] + (i % 2);
            Triplet_Add(massTriplet, global_i, global_i, localM[i][i]);  // Add only diagonal terms (lumped)
        }

        #ifdef LOC_MASS
        printf("Local lumped mass matrix for element %d:\n", e);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                printf("%12.5e ", localM[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        #endif
    }

    Triplet_Sort(massTriplet, 0, massTriplet->nz - 1);
    Triplet_Unique(massTriplet);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */    
void Mesh_Discretize_Dynamics(mesh *this, double *b, triplet *mass, triplet *p, triplet *dumping, double *dirichletValues, double (*psi_force)(double, double, double, double, double, double *, double *), double angle_degrees, const double eps1, const double eps2, int gravity) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    memset(b, 0, this->NPoints * sizeof(double) * 2); // Zero out the load vector

    assembleGlobalStiffnessMatrix(this, p, materials);
    computeConsistentMassMatrix(this, mass);
    Triplet_Damping(dumping, mass, p, eps1, eps2);

    if (gravity)
    {
        applyGravity(this, b, materials);
    }

    applyBoundaryConditions(this, b, dirichletValues, mass, psi_force, angle_degrees);
    applyBoundaryConditions(this, b, dirichletValues, p, psi_force, angle_degrees);
    applyBoundaryConditions(this, b, dirichletValues, dumping, psi_force, angle_degrees);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
//                                     OHYB NOSNIKU - 3D case
/* ------------------------------------------------------------------------------------------------- */
void calculateBMatrix3D(const gelement3D *K, double B[6][24], double G[3][8]) {
    int nbNodes = K->Nvertex;

    memset(B, 0, sizeof(double) * 6 * nbNodes * 3);

    for (int i = 0; i < nbNodes; i++) {
        B[0][3 * i + 0] = G[0][i]; // dN/dx
        B[1][3 * i + 1] = G[1][i]; // dN/dy
        B[2][3 * i + 2] = G[2][i]; // dN/dz

        B[3][3 * i + 0] = G[1][i]; // dN/dy
        B[3][3 * i + 1] = G[0][i]; // dN/dx

        B[4][3 * i + 1] = G[2][i]; // dN/dz
        B[4][3 * i + 2] = G[1][i]; // dN/dy

        B[5][3 * i + 0] = G[2][i]; // dN/dz
        B[5][3 * i + 2] = G[0][i]; // dN/dx
    }
}
/* ------------------------------------------------------------------------------------------------- */
//                                     ORTHOTROPIC MATERIAL
/* ------------------------------------------------------------------------------------------------- */
void computeRotationMatrix(double angleX, double angleY, double angleZ, double A[3][3]) {
    double cx = cos(angleX * M_PI / 180.0), sx = sin(angleX * M_PI / 180.0);
    double cy = cos(angleY * M_PI / 180.0), sy = sin(angleY * M_PI / 180.0);
    double cz = cos(angleZ * M_PI / 180.0), sz = sin(angleZ * M_PI / 180.0);

    double Rx[3][3] = {{1, 0, 0}, {0, cx, -sx}, {0, sx, cx}};
    double Ry[3][3] = {{cy, 0, sy}, {0, 1, 0}, {-sy, 0, cy}};
    double Rz[3][3] = {{cz, -sz, 0}, {sz, cz, 0}, {0, 0, 1}};

    // Matrix multiplication: R = Rz * Ry * Rx
    double temp[3][3];
    // temp = Rz * Ry
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            temp[i][j] = 0;
            for (int k = 0; k < 3; k++)
                temp[i][j] += Rz[i][k] * Ry[k][j];
        }
    // R = temp * Rx
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            A[i][j] = 0;
            for (int k = 0; k < 3; k++)
                A[i][j] += temp[i][k] * Rx[k][j];
        }
}
/* ------------------------------------------------------------------------------------------------- */
void buildVoigtTransformationMatrix(const double A[3][3], double T[6][6]) {
    memset(T, 0, sizeof(double) * 36);

    // T11 block
    T[0][0] = A[0][0]*A[0][0];
    T[0][1] = A[0][1]*A[0][1];
    T[0][2] = A[0][2]*A[0][2];
    T[1][0] = A[1][0]*A[1][0];
    T[1][1] = A[1][1]*A[1][1];
    T[1][2] = A[1][2]*A[1][2];
    T[2][0] = A[2][0]*A[2][0];
    T[2][1] = A[2][1]*A[2][1];
    T[2][2] = A[2][2]*A[2][2];

    // T12 block
    T[0][3] = 2.0 * A[0][0]*A[0][1];
    T[0][4] = 2.0 * A[0][0]*A[0][2];
    T[0][5] = 2.0 * A[0][1]*A[0][2];
    T[1][3] = 2.0 * A[1][0]*A[1][1];
    T[1][4] = 2.0 * A[1][0]*A[1][2];
    T[1][5] = 2.0 * A[1][1]*A[1][2];
    T[2][3] = 2.0 * A[2][0]*A[2][1];
    T[2][4] = 2.0 * A[2][0]*A[2][2];
    T[2][5] = 2.0 * A[2][1]*A[2][2];

    // T21 block
    T[3][0] = A[0][0]*A[1][0];
    T[3][1] = A[0][1]*A[1][1];
    T[3][2] = A[0][2]*A[1][2];
    T[4][0] = A[0][0]*A[2][0];
    T[4][1] = A[0][1]*A[2][1];
    T[4][2] = A[0][2]*A[2][2];
    T[5][0] = A[1][0]*A[2][0];
    T[5][1] = A[1][1]*A[2][1];
    T[5][2] = A[1][2]*A[2][2];

    // T22 block
    T[3][3] = A[0][0]*A[1][1] + A[0][1]*A[1][0];
    T[3][4] = A[0][0]*A[1][2] + A[0][2]*A[1][0];
    T[3][5] = A[0][1]*A[1][2] + A[0][2]*A[1][1];
    T[4][3] = A[0][0]*A[2][1] + A[0][1]*A[2][0];
    T[4][4] = A[0][0]*A[2][2] + A[0][2]*A[2][0];
    T[4][5] = A[0][1]*A[2][2] + A[0][2]*A[2][1];
    T[5][3] = A[1][0]*A[2][1] + A[1][1]*A[2][0];
    T[5][4] = A[1][0]*A[2][2] + A[1][2]*A[2][0];
    T[5][5] = A[1][1]*A[2][2] + A[1][2]*A[2][1];
}
/* ------------------------------------------------------------------------------------------------- */
void transformElasticityMatrix(const double D_local[6][6], const double T[6][6], double D_global[6][6]) {
    double temp[6][6] = {0};

    // temp = T * D_local
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++)
                temp[i][j] += T[i][k] * D_local[k][j];

    // D_global = temp * T^T
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            D_global[i][j] = 0.0;
            for (int k = 0; k < 6; k++)
                D_global[i][j] += temp[i][k] * T[j][k];
        }
}
/* ------------------------------------------------------------------------------------------------- */
int invert6x6(const double A[6][6], double A_inv[6][6]) {
    double temp[6][12];
    memset(temp, 0, sizeof(temp));

    // Create [A | I] augmented matrix
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            temp[i][j] = A[i][j];
            temp[i][j + 6] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan elimination
    for (int i = 0; i < 6; ++i) {
        // Pivot
        double pivot = temp[i][i];
        if (fabs(pivot) < 1e-12)
            return -1;  // Singular

        for (int j = 0; j < 12; ++j)
            temp[i][j] /= pivot;

        // Eliminate
        for (int k = 0; k < 6; ++k) {
            if (k == i) continue;
            double factor = temp[k][i];
            for (int j = 0; j < 12; ++j)
                temp[k][j] -= factor * temp[i][j];
        }
    }

    // Extract inverse
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            A_inv[i][j] = temp[i][j + 6];

    return 0;
}
/* ------------------------------------------------------------------------------------------------- */
void computeElasticityMatrixOrthotropic(MaterialProperties mat, double D[6][6], double angleX, double angleY, double angleZ) {
    
    double E1 = mat.E1, E2 = mat.E2, E3 = mat.E3;
    double nu12 = mat.nu12, nu13 = mat.nu13, nu23 = mat.nu23;
    double G12 = mat.G12, G13 = mat.G13, G23 = mat.G23;

    double nu21 = nu12 * (E2 / E1);
    double nu31 = nu13 * (E3 / E1);
    double nu32 = nu23 * (E3 / E2);

    memset(D, 0, sizeof(double) * 6 * 6);

    double Ddenom = ( 1.0 - nu12 * nu21 - nu23 * nu32 - nu13 * nu31 - 2.0 * nu21 * nu32 * nu13 ) / ( E1 * E2 * E3 );

    D[0][0] = (1 - nu23 * nu32) / (E2 * E3 * Ddenom);
    D[1][1] = (1 - nu13 * nu31) / (E1 * E3 * Ddenom);
    D[2][2] = (1 - nu12 * nu21) / (E1 * E2 * Ddenom);

    D[0][1] = (nu21 + nu23 * nu31) / (E2 * E3 * Ddenom);
    D[0][2] = (nu31 + nu21 * nu32) / (E2 * E3 * Ddenom);

    D[1][0] = (nu21 + nu23 * nu31) / (E2 * E3 * Ddenom);
    D[1][2] = (nu32 + nu12 * nu31) / (E1 * E3 * Ddenom);

    D[2][0] = (nu31 + nu21 * nu32) / (E2 * E3 * Ddenom);
    D[2][1] = (nu32 + nu12 * nu31) / (E1 * E3 * Ddenom);

    D[3][3] = G23;
    D[4][4] = G13;
    D[5][5] = G12;
    
    if (angleX != 0.0 || angleY != 0.0 || angleZ != 0.0) {
        // printf("Angel triggered\n");
        double A[3][3];
        double T[6][6], D_rot[6][6];
        computeRotationMatrix(angleX, angleY, angleZ, A);
        buildVoigtTransformationMatrix(A, T);
        transformElasticityMatrix(D, T, D_rot);
        memcpy(D, D_rot, sizeof(double) * 36);
    }
}
/* ------------------------------------------------------------------------------------------------- */
void computeElasticityMatrix3D(MaterialProperties mat, double D[6][6], double angleX, double angleY, double angleZ) {
    if (mat.model == ORTHOTROPIC) {
        computeElasticityMatrixOrthotropic(mat, D, angleX, angleY, angleZ);
        // printf("Orthotropic\n");
    } else {
        double E = mat.E;
        double nu = mat.nu;
        double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
        double mu = E / (2 * (1 + nu));

        memset(D, 0, sizeof(double) * 6 * 6);
        D[0][0] = D[1][1] = D[2][2] = lambda + 2 * mu;
        D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = lambda;
        D[3][3] = D[4][4] = D[5][5] = mu;
        // printf("Standard\n");
    }
}
/* ------------------------------------------------------------------------------------------------- */
void computeLocalStiffnessMatrix3D(mesh *this, int elementIndex, double localK[24][24], double angleX, double angleY, double angleZ, MaterialProperties *materials) {
    gelement3D K;
    GetElement3D(this, &K, elementIndex);

    double D[6][6];
    MaterialProperties mat = materials[elementIndex];
    computeElasticityMatrix3D(mat, D, angleX, angleY, angleZ);

    memset(localK, 0, sizeof(double) * 24 * 24);

    int nBaseNodes = GetNumBasisFunctions(K.type);
    int nGaussPoints = GetNumGaussPoints(K.type);

    for (int gp = 0; gp < nGaussPoints; gp++) {
        double xi, eta, zeta, weight;
        GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);

        double B[6][24] = {0};
        double G[3][8] = {0};
        double J[3][3], J_inv[3][3];
        double det_J = 0.0;

        // Compute gradients and Jacobian info
        computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, J, J_inv, &det_J, G, xi, eta, zeta);
        calculateBMatrix3D(&K, B, G);

        // Integrate contribution at Gauss point
        double dV = det_J * weight;
        for (int i = 0; i < nBaseNodes * 3; i++) {
            for (int j = 0; j < nBaseNodes * 3; j++) {
                for (int a = 0; a < 6; a++) {
                    for (int b = 0; b < 6; b++) {
                        localK[i][j] += ( B[a][i] * D[a][b] * B[b][j] ) * dV;
                    }
                }
            }
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void assembleGlobalStiffnessMatrix3D(mesh *this, triplet *p, double angleX, double angleY, double angleZ, MaterialProperties *materials) {
    double localK[24][24];
    for (int k = this->NTriangles + this->NQuads; k < this->NElements; k++) {
        gelement3D K;
        GetElement3D(this, &K, k);
        computeLocalStiffnessMatrix3D(this, k, localK, angleX, angleY, angleZ, materials);

        int nBaseNodes = GetNumBasisFunctions(K.type);

        for (int i = 0; i < nBaseNodes * 3; i++) {
            for (int j = 0; j < nBaseNodes * 3; j++) {
                int global_i = 3 * K.idxNode[i / 3] + (i % 3);
                int global_j = 3 * K.idxNode[j / 3] + (j % 3);
                Triplet_Add(p, global_i, global_j, localK[i][j]);
            }
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void applyGravity3D(mesh *this, double *b, MaterialProperties *materials) {
    double gravity = -9.81;
    for (int k = this->NTriangles + this->NQuads; k < this->NElements; k++) {
        gelement3D K;
        GetElement3D(this, &K, k);
        MaterialProperties mat = materials[k];

        int nGaussPoints = GetNumGaussPoints(K.type);
        double elementVolume = 0.0;

        for (int gp = 0; gp < nGaussPoints; gp++) {
            double xi, eta, zeta, weight;
            GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);

            double J[3][3], J_inv[3][3];
            double det_J = 0.0;
            double G[3][8] = {0};

            computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, J, J_inv, &det_J, G, xi, eta, zeta);

            elementVolume += det_J * weight;
        }

        double elementWeight = mat.rho * elementVolume * gravity;
        double loadPerNode = elementWeight / K.Nvertex;
        // printf("NVertex: %d\n",K.Nvertex);
        // printf("Rho %g \t gravity: %g \t K.detJ: %g\n", mat.rho, gravity, elementVolume);

        // printf("elementWeight = %g\n",elementWeight);

        for (int i = 0; i < K.Nvertex; i++) {
            b[3 * K.idxNode[i] + 1] += loadPerNode;
            // printf("Load per node: %g\n", loadPerNode);
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void applyBoundaryConditions3D(mesh *this, double *b, double *dirichletValues, triplet *p, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg) {
    
    for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
        if (this->ElementMark[elem] == NEUMANN_LOAD) {

            gelement3D face;
            GetElement3D(this, &face, elem);

            int nVertices = GetNumBasisFunctions(face.type);
            int nGaussPoints = GetNumGaussPoints(face.type);

            for (int gp = 0; gp < nGaussPoints; gp++) {
                double xi, eta, zeta, weight;
                GetGaussPoint(face.type, gp, &xi, &eta, &zeta, &weight);

                double gradN_xi[8] = {0.0}, gradN_eta[8] = {0.0};
                GetReferenceGradients(face.type, gradN_xi, gradN_eta, NULL, xi, eta, zeta);

                double J_face = 0.0, normal_gp[3] = {0.0};
                computeSurfaceJacobianAndNormal(face.Ver, nVertices, gradN_xi, gradN_eta, &J_face, normal_gp, NULL);

                double x, y, z;
                ComputePhysicalCoordinates(&x, &y, &z, &face, xi, eta, zeta);

                // Vector traction instead of scalar psi
                double traction[3] = {0.0};
                psi_force_3D(x, y, z, normal_gp[0], normal_gp[1], normal_gp[2], F_magnitude, azimuth_deg, elevation_deg, &traction[0], &traction[1], &traction[2]);

                double N[8];
                ComputeBasisFunctions(face.type, N, xi, eta, zeta, nVertices);

                // Assemble vector traction into global RHS
                for (int i = 0; i < nVertices; i++) {
                    int node = face.idxNode[i];
                    for (int d = 0; d < 3; d++) { // for u_x, u_y, u_z components
                        b[3 * node + d] += J_face * weight * N[i] * traction[d];
                        // printf("Tranction[%d] = %g\n", d, traction[d]);
                    }
                }
            }
        }
    }   
    
    for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
        if (this->ElementMark[elem] == DIRICHLET_BC) {
            
            gelement3D face;
            GetElement3D(this, &face, elem);

            int nVertices = GetNumBasisFunctions(face.type);

            for (int i = 0; i < nVertices; i++)
            {
                int node = face.idxNode[i];
                
                int ux = 3 * node;
                int uy = 3 * node + 1;
                int uz = 3 * node + 2; 

                for (int nz = 0; nz < p->nz; nz++) {
                    if (p->tripletI[nz] == ux || p->tripletI[nz] == uy || p->tripletI[nz] == uz || p->J[nz] == ux || p->J[nz] == uy || p->J[nz] == uz) {
                        p->VAL[nz] = 0.0;
                    }
                }
    
                Triplet_Add(p, ux, ux, 1.0);
                Triplet_Add(p, uy, uy, 1.0);
                Triplet_Add(p, uz, uz, 1.0);
    
                b[ux] = dirichletValues[ux];
                b[uy] = dirichletValues[uy];
                b[uz] = dirichletValues[uz];
            }
        }
    }

    // Artificially constrain all Z-DOFs to simulate pure XY bending
    // for (int i = 0; i < this->NPoints; i++) {
    //     int uz = 3 * i + 2; // Z-dof index
    //     Triplet_Add(p, uz, uz, 1.0);
    //     b[uz] = 0.0;
    // }


    Triplet_Sort(p, 0, p->nz - 1);
    Triplet_Unique(p);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Statics3D(mesh *this, double *b, triplet *p, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, int gravity, double angleX, double angleY, double angleZ) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    memset(b, 0, this->NPoints * sizeof(double) * 3);

    assembleGlobalStiffnessMatrix3D(this, p, angleX, angleY, angleZ, materials);

    if (gravity) {
        applyGravity3D(this, b, materials);
        // for (int i = 0; i < this->NPoints; i++)
        // {
        //     printf("b[%d] = %g\n", 3 * i + 1, b[3 * i + 1]);
        // }
    }

    applyBoundaryConditions3D(this, b, dirichletValues, p, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Statics3D_Orthotropic(mesh *this, double *b, triplet *p, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, int gravity, double angleX, double angleY, double angleZ) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    memset(b, 0, this->NPoints * sizeof(double) * 3);

    assembleGlobalStiffnessMatrix3D(this, p, angleX, angleY, angleZ, materials);
    
    // Apply gravity load (if specified)
    if (gravity) {
        applyGravity3D(this, b, materials);
    }

    // Apply boundary conditions (Dirichlet, Neumann, etc.)
    applyBoundaryConditions3D(this, b, dirichletValues, p, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);

    // Free the allocated memory for materials
    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
//                                     DYNAMICS - 3D case
/* ------------------------------------------------------------------------------------------------- */
void computeConsistentMassMatrix3D(mesh *this, triplet *massTriplet, MaterialProperties *materials) {
    for (int e = 0; e < this->NElements; e++) {
        gelement3D K;
        GetElement3D(this, &K, e);

        int nBaseNodes = GetNumBasisFunctions(K.type);
        double localM[24][24] = {0};
        double rho = materials[e].rho;
        int nGaussPoints = GetNumGaussPoints(K.type);

        for (int gp = 0; gp < nGaussPoints; gp++) {
            double xi, eta, zeta, weight;
            GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);
            double N[8] = {0};
            ComputeBasisFunctions(K.type, N, xi, eta, zeta, nBaseNodes);

            double J[3][3], J_inv[3][3], det_J;
            double G[3][8];
            computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, J, J_inv, &det_J, G, xi, eta, zeta);

            double dV = det_J * weight;

            for (int i = 0; i < nBaseNodes; i++) {
                for (int j = 0; j < nBaseNodes; j++) {
                    double m = rho * N[i] * N[j] * dV;
                    localM[3 * i + 0][3 * j + 0] += m;
                    localM[3 * i + 1][3 * j + 1] += m;
                    localM[3 * i + 2][3 * j + 2] += m;
                }
            }
        }

        for (int i = 0; i < nBaseNodes * 3; i++) {
            for (int j = 0; j < nBaseNodes * 3; j++) {
                int global_i = 3 * K.idxNode[i / 3] + (i % 3);
                int global_j = 3 * K.idxNode[j / 3] + (j % 3);
                Triplet_Add(massTriplet, global_i, global_j, localM[i][j]);
            }
        }
    }

    Triplet_Sort(massTriplet, 0, massTriplet->nz - 1);
    Triplet_Unique(massTriplet);
}
/* ------------------------------------------------------------------------------------------------- */
void assembleGlobalLumpedMassMatrix3D(mesh *this, triplet *massTriplet, MaterialProperties *materials) {
    for (int e = 0; e < this->NElements; e++) {
        gelement3D K;
        GetElement3D(this, &K, e);

        int nBaseNodes = GetNumBasisFunctions(K.type);
        double rho = materials[e].rho;
        int nGaussPoints = GetNumGaussPoints(K.type);

        double N[8];  // max number of nodes; adjust or malloc if needed
        for (int gp = 0; gp < nGaussPoints; gp++) {
            double xi, eta, zeta, weight;
            GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);

            ComputeBasisFunctions(K.type, N, xi, eta, zeta, nBaseNodes);

            double J[3][3], J_inv[3][3], det_J;
            double G[3][8];
            computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, J, J_inv, &det_J, G, xi, eta, zeta);

            double dV = det_J * weight;

            for (int i = 0; i < nBaseNodes; i++) {
                double lump_mass = rho * N[i] * dV;
                int node = K.idxNode[i];

                Triplet_Add(massTriplet, 3 * node + 0, 3 * node + 0, lump_mass);
                Triplet_Add(massTriplet, 3 * node + 1, 3 * node + 1, lump_mass);
                Triplet_Add(massTriplet, 3 * node + 2, 3 * node + 2, lump_mass);
            }
        }
    }

    Triplet_Sort(massTriplet, 0, massTriplet->nz - 1);
    Triplet_Unique(massTriplet);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Dynamics3D(mesh *this, double *b, triplet *mass, triplet *p, triplet *damping, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, const double eps1, const double eps2, int gravity, double angleX, double angleY, double angleZ) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    memset(b, 0, this->NPoints * sizeof(double) * 3);

    assembleGlobalStiffnessMatrix3D(this, p, angleX, angleY, angleZ, materials);
    computeConsistentMassMatrix3D(this, mass, materials);
    Triplet_Damping(damping, mass, p, eps1, eps2);

    if (gravity) {
        applyGravity3D(this, b, materials);
    }

    applyBoundaryConditions3D(this, b, dirichletValues, mass, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);
    applyBoundaryConditions3D(this, b, dirichletValues, p, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);
    applyBoundaryConditions3D(this, b, dirichletValues, damping, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Dynamics3D_Orthotropic(mesh *this, double *b, triplet *mass, triplet *p, triplet *damping, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg, const double eps1, const double eps2, int gravity, double angleX, double angleY, double angleZ) {
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    memset(b, 0, this->NPoints * sizeof(double) * 3);

    assembleGlobalStiffnessMatrix3D(this, p, angleX, angleY, angleZ, materials);
    computeConsistentMassMatrix3D(this, mass, materials);
    // assembleGlobalLumpedMassMatrix3D(this, mass, materials);
    Triplet_Damping(damping, mass, p, eps1, eps2);

    if (gravity) {
        applyGravity3D(this, b, materials);
    }

    applyBoundaryConditions3D(this, b, dirichletValues, mass, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);
    applyBoundaryConditions3D(this, b, dirichletValues, p, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);
    applyBoundaryConditions3D(this, b, dirichletValues, damping, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
//                                     NON-LINEAR ELASTICITY
/* ------------------------------------------------------------------------------------------------- */
void computeLameParameters(MaterialProperties mat, double *lambda, double *mu) {
    
    *lambda = (mat.E * mat.nu) / ((1 + mat.nu) * (1 - 2 * mat.nu));
    *mu = mat.E / (2 * (1 + mat.nu));

}
/* ------------------------------------------------------------------------------------------------- */
//                                      SVK-Specific
/* ------------------------------------------------------------------------------------------------- */
void computeSVKStress(MaterialProperties mat, double F[3][3], double sigma[3][3]) {
    double lambda, mu;
    computeLameParameters(mat, &lambda, &mu);

    // Step 1: Compute right Cauchy-Green deformation tensor: C = Fᵀ * F
    double C[3][3] = {{0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                C[i][j] += F[k][i] * F[k][j];  // Note: Fᵀ * F
            }
        }
    }

    // Step 2: Compute Green-Lagrange strain tensor: E = 0.5 * (C - I)
    double E[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            E[i][j] = 0.5 * (C[i][j] - (i == j ? 1.0 : 0.0));
        }
    }

    // Step 3: Compute 2nd Piola-Kirchhoff stress: S = λ * tr(E) * I + 2μ * E
    double trace_E = E[0][0] + E[1][1] + E[2][2];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sigma[i][j] = lambda * trace_E * (i == j ? 1.0 : 0.0) + 2.0 * mu * E[i][j];
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void computeTangentModulusSVK(MaterialProperties mat, double C[6][6]) {
    double lambda, mu;
    computeLameParameters(mat, &lambda, &mu);
    memset(C, 0, sizeof(double) * 6 * 6);

    C[0][0] = C[1][1] = C[2][2] = lambda + 2.0 * mu;
    C[0][1] = C[0][2] = C[1][0] = C[1][2] = C[2][0] = C[2][1] = lambda;
    C[3][3] = C[4][4] = C[5][5] = mu;
}
/* ------------------------------------------------------------------------------------------------- */
//                                      NEOHOOKE-Specific
/* ------------------------------------------------------------------------------------------------- */
int voigt_index(int i, int j) {
    if (i == j) return i; // xx, yy, zz → 0,1,2
    if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3; // yz
    if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4; // zx
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 5; // xy
    return -1; // error
}
/* ------------------------------------------------------------------------------------------------- */
void computeNeoHookeanStress(MaterialProperties mat, double F[3][3], double sigma[3][3]) {
    double lambda, mu;
    computeLameParameters(mat, &lambda, &mu);

    if (mat.nu >= 0.4999) {
        fprintf(stderr, "[Error] Neo-Hookean model not valid for near-incompressible material (ν = %.5f). Use a mixed formulation or compressible model.\n", mat.nu);
        exit(EXIT_FAILURE);
    }

    // Compute J = det(F)
    double J = F[0][0]*(F[1][1]*F[2][2] - F[1][2]*F[2][1])
             - F[0][1]*(F[1][0]*F[2][2] - F[1][2]*F[2][0])
             + F[0][2]*(F[1][0]*F[2][1] - F[1][1]*F[2][0]);

    if (J <= 0.0) {
        fprintf(stderr, "[Error] Invalid deformation gradient: det(F) = %e <= 0\n", J);
        exit(EXIT_FAILURE);
    }

    // Compute C = Fᵀ * F
    double C[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C[i][j] += F[k][i] * F[k][j];

    // Compute C⁻¹
    double detC = C[0][0]*(C[1][1]*C[2][2] - C[1][2]*C[2][1])
                - C[0][1]*(C[1][0]*C[2][2] - C[1][2]*C[2][0])
                + C[0][2]*(C[1][0]*C[2][1] - C[1][1]*C[2][0]);

    if (detC <= 1e-12) {
        fprintf(stderr, "[Error] det(C) too small: %e. Aborting to avoid singular inversion.\n", detC);
        exit(EXIT_FAILURE);
    }                

    double Cinv[3][3] = {0};
    double invDet = 1.0 / detC;

    Cinv[0][0] =  (C[1][1]*C[2][2] - C[1][2]*C[2][1]) * invDet;
    Cinv[0][1] = -(C[0][1]*C[2][2] - C[0][2]*C[2][1]) * invDet;
    Cinv[0][2] =  (C[0][1]*C[1][2] - C[0][2]*C[1][1]) * invDet;
    Cinv[1][0] = -(C[1][0]*C[2][2] - C[1][2]*C[2][0]) * invDet;
    Cinv[1][1] =  (C[0][0]*C[2][2] - C[0][2]*C[2][0]) * invDet;
    Cinv[1][2] = -(C[0][0]*C[1][2] - C[0][2]*C[1][0]) * invDet;
    Cinv[2][0] =  (C[1][0]*C[2][1] - C[1][1]*C[2][0]) * invDet;
    Cinv[2][1] = -(C[0][0]*C[2][1] - C[0][1]*C[2][0]) * invDet;
    Cinv[2][2] =  (C[0][0]*C[1][1] - C[0][1]*C[1][0]) * invDet;

    // Symmetrize
    for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 3; j++)
            Cinv[i][j] = Cinv[j][i] = 0.5 * (Cinv[i][j] + Cinv[j][i]);

    // Compute 2nd Piola–Kirchhoff stress
    double lnJ = log(J);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            sigma[i][j] = mu * ((i == j ? 1.0 : 0.0) - Cinv[i][j]) + lambda * lnJ * Cinv[i][j];
}
/* ------------------------------------------------------------------------------------------------- */
void computeTangentModulusNeoHooke(const double F[3][3], MaterialProperties mat, double C_voigt[6][6]) {
    double lambda, mu;
    computeLameParameters(mat, &lambda, &mu);

    // Step 1: Compute right Cauchy-Green tensor C = Fᵀ * F
    double C[3][3] = {0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C[i][j] += F[k][i] * F[k][j];

    // Step 2: Compute determinant of C = J^2
    double detC = C[0][0]*(C[1][1]*C[2][2] - C[1][2]*C[2][1])
                - C[0][1]*(C[1][0]*C[2][2] - C[1][2]*C[2][0])
                + C[0][2]*(C[1][0]*C[2][1] - C[1][1]*C[2][0]);

    if (detC <= 0.0) {
        printf("Warning: det(C) = %e <= 0, cannot compute exact tangent.\n", detC);
        memset(C_voigt, 0, sizeof(double) * 36);
        return;
    }

    // double logJ = 0.5 * log(detC);

    // Step 3: Compute inverse of C
    double C_inv[3][3];
    double invDet = 1.0 / detC;

    C_inv[0][0] =  (C[1][1]*C[2][2] - C[1][2]*C[2][1]) * invDet;
    C_inv[0][1] = -(C[0][1]*C[2][2] - C[0][2]*C[2][1]) * invDet;
    C_inv[0][2] =  (C[0][1]*C[1][2] - C[0][2]*C[1][1]) * invDet;
    C_inv[1][0] = -(C[1][0]*C[2][2] - C[1][2]*C[2][0]) * invDet;
    C_inv[1][1] =  (C[0][0]*C[2][2] - C[0][2]*C[2][0]) * invDet;
    C_inv[1][2] = -(C[0][0]*C[1][2] - C[0][2]*C[1][0]) * invDet;
    C_inv[2][0] =  (C[1][0]*C[2][1] - C[1][1]*C[2][0]) * invDet;
    C_inv[2][1] = -(C[0][0]*C[2][1] - C[0][1]*C[2][0]) * invDet;
    C_inv[2][2] =  (C[0][0]*C[1][1] - C[0][1]*C[1][0]) * invDet;

    // Step 4: Symmetrize C-1
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            C_inv[i][j] = C_inv[j][i] = 0.5 * (C_inv[i][j] + C_inv[j][i]);

    // Step 5: Build tangent modulus in Voigt notation
    memset(C_voigt, 0, sizeof(double) * 6 * 6);

    for (int I = 0; I < 3; ++I) {
        for (int J = 0; J < 3; ++J) {
            int a = voigt_index(I, J);
            if (a < 0) continue;

            for (int K = 0; K < 3; ++K) {
                for (int L = 0; L < 3; ++L) {
                    int b = voigt_index(K, L);
                    if (b < 0) continue;

                    double Cijkl = lambda * C_inv[I][J] * C_inv[K][L]
                                 + (mu - lambda * log(sqrt(detC))) * (C_inv[I][K] * C_inv[J][L] + C_inv[I][L] * C_inv[J][K]);

                    C_voigt[a][b] += Cijkl;
                }
            }
        }
    }

    for(int a = 0; a < 6; ++a) {
        double fa = (a>=3) ? 2.0 : 1.0;        /* shear row scale */
        for(int b = 0; b < 6; ++b) {
            double fb = (b>=3) ? 2.0 : 1.0;    /* shear col scale */
            C_voigt[a][b] *= fa * fb;
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
//                                      SHARED FUNCTIONS
/* ------------------------------------------------------------------------------------------------- */
void computeStressNonLinear(MaterialProperties mat, double F[3][3], double sigma[3][3]) {
    if (mat.model == SVK) {
        computeSVKStress(mat, F, sigma);
        // printf("Saint-Venant Kirchhoff (SVK) model used\n");
    } else if (mat.model == NEOHOOKE) {
        computeNeoHookeanStress(mat, F, sigma);
        // printf("NeoHooke model used\n");
    }
}
/* ------------------------------------------------------------------------------------------------- */
void computeDeformationGradient(const gelement3D *K, const double *u_global, double F[3][3], double G[3][8]) {
    // Identity initialization
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // F = I + grad_U
    for (int a = 0; a < K->Nvertex; a++) {
        int node = K->idxNode[a];

        double ua[3] = {
            u_global[3 * node + 0], //u_x
            u_global[3 * node + 1], //u_y
            u_global[3 * node + 2]  //u_z
        };

        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++)
            {
                F[i][j] += ua[i] * G[j][a];
            }
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void ApplyDirichletBCsFromFaces3D(mesh *M, triplet *stiff, double *residual, double *f_ext, double *P_int, double *u_global, double *dirichletValues) {
    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D face;
            GetElement3D(M, &face, elem);

            int nVertices = GetNumBasisFunctions(face.type);
            for (int i = 0; i < nVertices; i++) {
                int node = face.idxNode[i];
                for (int d = 0; d < 3; d++) {
                    int dof = 3 * node + d;
                    for (int nz = 0; nz < stiff->nz; nz++) {
                        if (stiff->tripletI[nz] == dof || stiff->J[nz] == dof) {
                            stiff->VAL[nz] = 0.0;
                        }
                    }
                    Triplet_Add(stiff, dof, dof, 1.0);
                    residual[dof] = u_global[dof] - dirichletValues[dof];
                    f_ext[dof] = 0.0;
                    P_int[dof] = 0.0;

                }
            }
        }
    }

    Triplet_Sort(stiff, 0, stiff->nz - 1);
    Triplet_Unique(stiff);
}
/* ------------------------------------------------------------------------------------------------- */
void ApplyDirichletBCsToMassMatrix(mesh *M, triplet *mass) {
    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D face;
            GetElement3D(M, &face, elem);

            int nVertices = GetNumBasisFunctions(face.type);
            for (int i = 0; i < nVertices; i++) {
                int node = face.idxNode[i];
                for (int d = 0; d < 3; d++) {
                    int dof = 3 * node + d;
                    for (int nz = 0; nz < mass->nz; nz++) {
                        if (mass->tripletI[nz] == dof || mass->J[nz] == dof) {
                            mass->VAL[nz] = 0.0;
                        }
                    }
                    Triplet_Add(mass, dof, dof, 1.0);  // Replace with identity entry
                }
            }
        }
    }

    Triplet_Sort(mass, 0, mass->nz - 1);
    Triplet_Unique(mass);
}
/* ------------------------------------------------------------------------------------------------- */
void computeInternalLocalResidualForceNonlinear3D(mesh *this, int elementIndex, double *u_global, double *r_local, MaterialProperties *mat) {
    gelement3D K;
    GetElement3D(this, &K, elementIndex);

    int nBaseNodes = GetNumBasisFunctions(K.type);
    int ndof = 3 * nBaseNodes;
    int nGauss = GetNumGaussPoints(K.type);

    memset(r_local, 0, sizeof(double) * ndof);

    for (int gp = 0; gp < nGauss; gp++) {
        double xi, eta, zeta, weight;
        GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);

        double B[6][24] = {{0}}; // MAX_DOF = 3 * MAX_NODES, defined globally
        double G[3][8] = {{0}};
        double J[3][3], J_inv[3][3], det_J;
        double sigma[3][3], F[3][3];

        computeJacobianAndInverse(K.Ver, nBaseNodes, K.type, J, J_inv, &det_J, G, xi, eta, zeta);
        computeDeformationGradient(&K, u_global, F, G);
        computeStressNonLinear(*mat, F, sigma);
        calculateBMatrix3D(&K, B, G);

        double dV = det_J * weight;
        if (mat->model == NEOHOOKE || mat->model == SVK) 
        {
            double sigma_voigt[6] = {
                sigma[0][0],                 /* σxx */
                sigma[1][1],                 /* σyy */
                sigma[2][2],                 /* σzz */
                sigma[0][1],                 /* τxy = σxy */
                sigma[1][2],                 /* τyz = σyz */
                sigma[2][0]                  /* τzx = σzx */ // 2.0 * 
            };
            for (int i = 0; i < ndof; i++) {
                for (int a = 0; a < 6; a++) {
                    r_local[i] += B[a][i] * sigma_voigt[a] * dV;
                }
            }
        } else {
            printf("Material model is nethier NeoHookean or SVK\n");
            exit(EXIT_FAILURE);
        } 
    }
}
/* ------------------------------------------------------------------------------------------------- */
void assembleInternalGlobalResidualForceNonlinear3D(mesh *M, double *residual, double *u_global, MaterialProperties *materials) {
    memset(residual, 0, sizeof(double) * 3 * M->NPoints);

    for (int e = M->NTriangles + M->NQuads; e < M->NElements; e++) {
        gelement3D K;
        GetElement3D(M, &K, e);
        int nBaseNodes = GetNumBasisFunctions(K.type);
        int ndof = 3 * nBaseNodes;
        double r_local[24];

        computeInternalLocalResidualForceNonlinear3D(M, e, u_global, r_local, &materials[e]);

        for (int i = 0; i < ndof; i++) {
            int global_i = 3 * K.idxNode[i / 3] + (i % 3);
            residual[global_i] += r_local[i];
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void computeExternalResidualForceNonlinear3D(mesh *M, double *f_ext, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_magnitude, double azimuth_deg, double elevation_deg) {
    memset(f_ext, 0, sizeof(double) * 3 * M->NPoints);
    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == NEUMANN_LOAD) {
            gelement3D face;
            GetElement3D(M, &face, elem);

            int nVertices = GetNumBasisFunctions(face.type);
            int nGauss = GetNumGaussPoints(face.type);

            for (int gp = 0; gp < nGauss; gp++) {
                double xi, eta, zeta, weight;
                GetGaussPoint(face.type, gp, &xi, &eta, &zeta, &weight);

                double gradN_xi[8], gradN_eta[8];
                GetReferenceGradients(face.type, gradN_xi, gradN_eta, NULL, xi, eta, zeta);

                double J_face, normal[3];
                computeSurfaceJacobianAndNormal(face.Ver, nVertices, gradN_xi, gradN_eta, &J_face, normal, NULL);

                double x, y, z;
                ComputePhysicalCoordinates(&x, &y, &z, &face, xi, eta, zeta);

                double traction[3] = {0.0};
                psi_force_3D(x, y, z, normal[0], normal[1], normal[2], F_magnitude, azimuth_deg, elevation_deg, &traction[0], &traction[1], &traction[2]);

                double N[8];
                ComputeBasisFunctions(face.type, N, xi, eta, zeta, nVertices);

                for (int i = 0; i < nVertices; i++) {
                    int node = face.idxNode[i];
                    for (int d = 0; d < 3; d++) {
                        f_ext[3 * node + d] += J_face * weight * N[i] * traction[d];
                    }
                }
            }
        }
    }
    // int neumann_faces = 0;
    // for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
    //     if (M->ElementMark[elem] == NEUMANN_LOAD)
    //         ++neumann_faces;
    // }
    // printf("‣ Found %d Neumann boundary faces.\n", neumann_faces);
}
/* ------------------------------------------------------------------------------------------------- */
void computeLocalTangentStiffnessMatrixNonlinear3D(mesh *this, int elementIndex, double *u_global, double localK[24][24], MaterialProperties *materials) {
    gelement3D K;
    GetElement3D(this, &K, elementIndex);

    MaterialProperties mat = materials[elementIndex];

    memset(localK, 0, sizeof(double) * 24 * 24);

    int nBaseNodes = GetNumBasisFunctions(K.type);
    int nGauss = GetNumGaussPoints(K.type);

    for (int gp = 0; gp < nGauss; gp++) {
        double xi, eta, zeta, weight;
        GetGaussPoint(K.type, gp, &xi, &eta, &zeta, &weight);
        
        double B[6][24] = {0};
        double G[3][8] = {0};
        double J[3][3], Jinv[3][3];
        double det_J = 0.0;
        double F[3][3];
        double C[6][6];
        
        computeJacobianAndInverse(K.Ver, K.Nvertex, K.type, J, Jinv, &det_J, G, xi, eta, zeta);
        computeDeformationGradient(&K, u_global, F, G);

        double S[3][3];                     /* 2nd PK stress*/
        computeStressNonLinear(mat, F, S);  /* fill S = σ */

        if (mat.model == SVK) {
            computeTangentModulusSVK(mat, C);
            // printf("SVK triggered\n");
        } else if (mat.model == NEOHOOKE) {
            computeTangentModulusNeoHooke(F, mat, C);
        }

        calculateBMatrix3D(&K, B, G);

        double dV = det_J * weight;
        for (int i = 0; i < nBaseNodes * 3; i++) {
            for (int j = 0; j < nBaseNodes * 3; j++) {
                for (int a = 0; a < 6; a++) {
                    for (int b = 0; b < 6; b++) {
                        localK[i][j] += B[a][i] * C[a][b] * B[b][j] * dV;
                    }
                }
            }
        }

        /* Gσ matrix :   Gσ_{ia} = n_j σ_{ji} φ_{,a}^j       (Hughes, eq. 20.3.37) */
        double Gsig[3][24] = {{0}};
        for (int a = 0; a < nBaseNodes; ++a)
        {
            for (int i = 0; i < 3; ++i)
            {
                /* σ_{ji} * φ_,j */
                for (int j = 0; j < 3; ++j)
                {
                    Gsig[i][3 * a + j] += S[j][i] * G[j][a];
                }                
            }
        }

        for (int i = 0; i < 3 * nBaseNodes; ++i) {
            for (int j = 0; j < 3 * nBaseNodes; ++j) {
                for (int k = 0; k < 3;  ++k) {         /* σ acts only on k∈{x,y,z} */
                    localK[i][j] += Gsig[k][i] * G[k][j / 3] * ((j % 3) == k) * dV;
                }
            }
        }
    }    
}
/* ------------------------------------------------------------------------------------------------- */
void assembleGlobalTangentStiffnessMatrixNonlinear3D(mesh *this, triplet *stiff, double *u_global, MaterialProperties *materials) {
    double localK[24][24];
    for (int k = this->NTriangles + this->NQuads; k < this->NElements; k++) {
        gelement3D K;
        GetElement3D(this, &K, k);
        computeLocalTangentStiffnessMatrixNonlinear3D(this, k, u_global, localK, materials);

        int nBaseNodes = GetNumBasisFunctions(K.type);

        for (int i = 0; i < nBaseNodes * 3; i++) {
            for (int j = 0; j < nBaseNodes * 3; j++) {
                int global_i = 3 * K.idxNode[i / 3] + (i % 3);
                int global_j = 3 * K.idxNode[j / 3] + (j % 3);
                Triplet_Add(stiff, global_i, global_j, localK[i][j]);
            }
        }
    }
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_NonlinearStatics3D(mesh *M, triplet *stiff, double *residual, double *f_ext, double *P_int, double *u_global, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg)
{
    MaterialProperties *materials = malloc(M->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(M, materials);

    // Assemble residual
    assembleInternalGlobalResidualForceNonlinear3D(M, P_int, u_global, materials); // P(u)

    // Assemble tangent stiffness
    assembleGlobalTangentStiffnessMatrixNonlinear3D(M, stiff, u_global, materials); //Kt

    // Compute external forces
    // double *f_ext = calloc(3 * M->NPoints, sizeof(double));
    computeExternalResidualForceNonlinear3D(M, f_ext, psi_force_3D, F_mag, azimuth_deg, elev_deg); // f

    // double f_sum = 0.0;
    // for (int i = 0; i < 3 * M->NPoints; ++i)
    //     f_sum += fabs(f_ext[i]);
    // printf("‣ Total applied external force = %.12e\n", f_sum);

    for (int i = 0; i < 3 * M->NPoints; i++)
        residual[i] = -f_ext[i] + P_int[i]; // f - P(u)

    // Apply Dirichlet BCs
    ApplyDirichletBCsFromFaces3D(M, stiff, residual, f_ext, P_int, u_global, dirichletValues);

    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_Mass_Matrix(mesh *this, triplet *mass)
{
    MaterialProperties *materials = malloc(this->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(this, materials);

    computeConsistentMassMatrix3D(this, mass, materials);
    ApplyDirichletBCsToMassMatrix(this, mass);
    
    free(materials);
}
/* ------------------------------------------------------------------------------------------------- */
void Mesh_Discretize_NonlinearDynamics3D(mesh *M, triplet *stiff, double *residual, double *f_ext, double *P_int, double *u_global, double *dirichletValues, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg)
{
    MaterialProperties *materials = malloc(M->NElements * sizeof(MaterialProperties));
    initializeMaterialProperties(M, materials);

    // Internal forces (P_int)
    assembleInternalGlobalResidualForceNonlinear3D(M, P_int, u_global, materials); // P(u)

    // Tangent stiffness matrix (K_t)
    assembleGlobalTangentStiffnessMatrixNonlinear3D(M, stiff, u_global, materials);

    computeExternalResidualForceNonlinear3D(M, f_ext, psi_force_3D, F_mag, azimuth_deg, elev_deg); // f

    for (int i = 0; i < 3 * M->NPoints; i++)
        residual[i] = -f_ext[i] + P_int[i]; // f - P(u)

    ApplyDirichletBCsFromFaces3D(M, stiff, P_int, residual, P_int, u_global, dirichletValues);

    free(materials);
}
