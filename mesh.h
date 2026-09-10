/* ----------------------------------------------------------------------------------------- */
#ifndef _MESH_H
#define _MESH_H
/* ----------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------- */
#define INSIDE 0
#define DIRICHLET_BC 1 
#define NEUMANN_BC 11
#define NEUMANN_LOAD 10
#define NEUMANN_FREE 11
#define ROBIN_BC 21
/* ----------------------------------------------------------------------------------------- */
#define STEEL 100
#define ALUMINUM 110
#define PVC 120
#define SILICONE 130
#define MAHAGONY 140
/* ----------------------------------------------------------------------------------------- */
#define pi 3.14159265358979323846
/* ----------------------------------------------------------------------------------------- */
#define g -9.81
/* ----------------------------------------------------------------------------------------- */
#define alfa 1.0
/* ----------------------------------------------------------------------------------------- */
#define BETA 0.25
#define GAMMA 0.5 //NEWMARK_
/* ----------------------------------------------------------------------------------------- */
typedef enum {
    INTEGRATION_VERTICES,
    INTEGRATION_MIDPOINTS,
    INTEGRATION_TEZISTE,
    INTEGRATION_GAUSSIAN
} IntegrationType;
/* ----------------------------------------------------------------------------------------- */
typedef enum { 
    SVK,
    NEOHOOKE,
    ORTHOTROPIC,
    LINEAR
} MaterialModel;

typedef struct {
    double E;   // Youngs modulus
    double nu;  // Poissons ratio
    double rho;  // Density
    MaterialModel model;

    //Orthotropic materials
    double E1, E2, E3;  
    double nu12, nu13, nu23;
    double G12, G13, G23;
    double angleX, angleY, angleZ;  // Material orientation in global coordinates (Euler angles)
} MaterialProperties;

/* ----------------------------------------------------------------------------------------- */
typedef struct {
    /* vertices */
    int     NPoints;
    double  *x, *y, *z;
    int     dim;

    int      NElements;               /**< Number of All 2D Elements */
    int      NTriangles,   NQuads;    /**< Number of Triangles & Quads */
    int      NTetrahedra;          /**< Number of tetrahedrons (3D) */
    int      NHexahedra;           /**< Number of hexahedrons (3D) */
    int      NWedges;              /**< Number of wedges/prisms (3D) */
    int      NPyramids;             /**< Number of pyramids (3D) */
    int      NBEdges;                 /**< Number of Boundary Edges */

    int     *ElementListOfVertices;   /**< List of elements vertices, ordered Triangles, Quads */
    short   *ElementMark;                  /**< Element mark */

    int     *EdgeVerA, *EdgeVerB;
    short   *EdgeMark;

    /* Zde treba zjistit z informaci nahore nasledujici: */
    short   *PointMark;                    /**< Point marks */
    short   *PointMark3D;                  /**< Point marks 3D */
    short   *NodeDegree;                   /**<  number of elements with the specified vertex */

    int    startQuads;                      /**<  starting point of Quads eleemtns */
    int    startTetra;
    int    startHex;
    int    startPyramid;
    int    startWedge;
    
    int     *PatchW;                       /**<  location of patch for the specified patch    */
    int     *PatchListOfElements;          /**<  list of elements for vertices in the patch   */

} mesh;                                    /**<  mesh structure */
/* ----------------------------------------------------------------------------------------- */
/* PUBLIC:                                                                                   */
/* ----------------------------------------------------------------------------------------- */
void   Mesh_Load(mesh *this, const char *fn);
void   Mesh_Free(mesh *this);
void   Mesh_Allocate(mesh *this);
double Mesh_GetDomainVolume(mesh *this); 
void   Mesh_ImportMSH(mesh *this, const char *fname);

void   PointMark(mesh *this);
void   Mesh_Interpolate(mesh *this, double f(double, double), double *uh);
void   IdentifyPointMark1(mesh *this, int **isDirichlet);

void exportDisplacementToVTK(const mesh *M, const double *displacement, const char *filename);
void exportNumericalSolutionToVTK(const mesh *M, const double *numericalSolution, const char *filename);
void saveSolutionToFile(const char* filename, double* solution, int N);
double* loadSolution(const char* filename, int* size);

double computeL2Error(mesh *M, double *numericalSolution, double (*analyticalSolution)(double, double));
double computeH1SeminormError(mesh *M, double *numericalSolution, double (*AnalyticalSolution)(double, double), void (*CalculateAnalyticalGradient)(double, double, double*, double*));
double computeL2Norm(double* solution_coarse, double* solution_fine_interpolated, int N);
double computeH1Seminorm(double* grad_solution_coarse, double* grad_solution_fine_interpolated, int N);

double computeL2Error3D(mesh *M, double *numericalSolution, double (*analyticalSolution)(double, double, double));
double computeH1SeminormError3D(mesh *M, double *numericalSolution, double (*AnalyticalSolution)(double, double, double), void (*CalculateAnalyticalGradient)(double, double, double, double*, double*, double*));

void initializeMaterialProperties(mesh *M, MaterialProperties* materials);

void ComputeBasisFunctions(int elementType, double *N, double xi, double eta, double zeta, int nNodes);
void GetReferenceGradients(int elementType, double *dN_xi, double *dN_eta, double *dN_zeta, double xi, double eta, double zeta);
int GetNumBasisFunctions(int elementType);

void computeResidualForceVector(mesh *this, double *displacement, double *residual, MaterialProperties *materials);


#endif