#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> //fabs, sin...
#include <assert.h> //assert
#include <cholmod.h>   // For CHOLMOD library
#include <umfpack.h>   // For UMFPACK library
#include <cs.h>        // For CSparse library
#include <petscksp.h> // For PETSc library

#include "mesh.h"
#include "constants.h"
#include "gelement.h"
#include "gelement3D.h"
#include "triplet.h"
#include "csr.h"
#include "integral.h"
#include "Mesh_GetBndrSide.h"
#include "saving.h"
#include "Poisson.h"
#include "dynamics.h"
#include "problems.h"

#include "Triplet_1.h"
#include "sparse_1.h"
#include "vector_1.h"

// #define POISSON2D
// #define STATICS2D
// #define DYNAMICS2D
// #define POISSON3D
// #define STATICS3D
// #define STATICS3D_ORTHO
// #define DYNAMICS3D
// #define DYNAMICS3D_ORTHO
#define NONLINEAR
// #define NONLINEAR_DYNAMICS

//------------------------------------------------------------------------------------------------------------
void Element_vertex_print(mesh *M, int Vertices);
void Mesh_Node_Coordinates_Print(mesh *M);
void Mesh_Classify_Nodes(mesh *M);
//------------------------------------------------------------------------------------------------------------

double f(double x, double y)
{
    return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

double f3D(double x, double y, double z)
{
    return 3 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) * 1.0 + 0.0;
}

double psi(double x, double y, double nx, double ny )
{
    return   nx * M_PI * cos(M_PI * x) * sin(M_PI * y) + ny * M_PI * sin(M_PI * x) * cos(M_PI * y);
}

double psi3D(double x, double y, double z, double nx, double ny, double nz )
{
    return nx * M_PI * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) +
           ny * M_PI * sin(M_PI * x) * cos(M_PI * y) * sin(M_PI * z) +
           nz * M_PI * sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * z);

    // return  M_PI * cos (M_PI * x) * sin (M_PI * y) * sin (M_PI * z) * 1.0 + 0.0;
}

double psi_force(double x, double y, double nx, double ny, double angle_degrees, double *force_x, double *force_y) {
    double angle_rad = angle_degrees * (M_PI / 180.0);
    double force_magnitude = 1e8;  // Assume constant force or any logic you want here

    // Compute force components in x and y directions based on angle and normals
    *force_x = force_magnitude * (cos(angle_rad) * nx - sin(angle_rad) * ny);
    *force_y = force_magnitude * (sin(angle_rad) * nx + cos(angle_rad) * ny);

    return sqrt((*force_x) * (*force_x) + (*force_y) * (*force_y));
}

double psi_force_3D(double x, double y, double z, double nx, double ny, double nz, double F_magnitude, double azimuth_deg, double elevation_deg, double *force_x, double *force_y, double *force_z) {
    double azimuth_rad = azimuth_deg * (M_PI / 180.0);
    double elevation_rad = elevation_deg * (M_PI / 180.0);

    double Fx = F_magnitude * cos(elevation_rad) * cos(azimuth_rad);
    double Fy = F_magnitude * cos(elevation_rad) * sin(azimuth_rad);
    double Fz = F_magnitude * sin(elevation_rad);

    // printf("Fx = %g \t Fy = %g \t Fz = %g \n", Fx, Fy, Fz);

    *force_x = Fx;
    *force_y = Fy;
    *force_z = Fz;

    return sqrt(Fx * Fx + Fy * Fy + Fz * Fz); // Just returning magnitude
}


double ur(double x, double y, double nx, double ny)
{
    // alfa - value defined in mesh.h
    return ( ( (nx * M_PI * cos(M_PI * x) * sin(M_PI * y) + ny * M_PI * sin(M_PI * x) * cos(M_PI * y))  / alfa ) + sin(M_PI * x) * sin(M_PI * y) );
}

double ur3D(double x, double y, double z, double nx, double ny, double nz)
{
    // alfa - value defined in mesh.h
    return ( ( (nx * M_PI * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) + ny * M_PI * sin(M_PI * x) * cos(M_PI * y) * sin(M_PI * z) + nz * M_PI * sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * z))  / alfa ) + sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) );
    // return M_PI * cos (M_PI * x) * sin (M_PI * y) * sin (M_PI * z) * 0.0 + 1.0;
}

//--------------------------------------------------------------------------------------------------------------
int main() {
#ifdef POISSON2D
    mesh M;
    triplet stiff;
    csr A_csr;
    double *b = NULL, *dirichletValues = NULL, *x = NULL;

    // Control parameters
    IntegrationType integrationType = INTEGRATION_VERTICES;
    int solverType = 1;    // 1 = UMFPACK, 2 = Iterative Solver
    double tolerance = 1e-9;

    Mesh_Load(&M, "Poisson_Geo/squareDIR.msh");

    InitializePoissonProblem(&M, &stiff, &A_csr, &b, &dirichletValues, &x, integrationType);
    ExecutePoissonAnalysis(&M, &A_csr, b, x, M.NPoints, solverType, tolerance);
    Poisson_ComputeNorms(&M, x);
    CleanupPoissonProblem(&stiff, &A_csr, b, dirichletValues, x);
#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef STATICS2D
    mesh M;
    triplet stiff;
    csr K_csr;
    double *b = NULL, *dirichletValues = NULL, *x = NULL;

    // Control parameters
    double angle_degrees = 180.0;
    int gravity = FALSE;
    int solverType = 1;    // 1 = UMFPACK, 2 = Iterative Solver
    double tolerance = 1e-9;

    Mesh_Load(&M, "Elasticity_Geo/MESH_2D/Beam2D_bend.msh");
    int nDOF = 2 * M.NPoints;

    InitializeStaticsProblem(&M, &stiff, &K_csr, &b, &dirichletValues, &x, angle_degrees, gravity);
    ExecuteStaticsAnalysis(&M, &K_csr, b, x, nDOF, solverType, tolerance);

    saveSparseMatrixToFile(&K_csr, "BEAM_solution/Stiffness_2D.dat");
    saveVectorToFile(b, 2 * M.NPoints, 1, "BEAM_solution/Load_vector_2D.txt");
    saveVectorToFile(x, 2 * M.NPoints, 1, "BEAM_solution/Solution_vector_2D.txt");

    saveVectorToFile(b, 2 * M.NPoints, 2, "BEAM_solution/Load_vector_2D_visual.txt");
    saveVectorToFile(x, 2 * M.NPoints, 2, "BEAM_solution/Solution_vector_2D_visual.txt");

    CleanupStaticsProblem(&stiff, &K_csr, b, dirichletValues, x);

#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef DYNAMICS2D
    mesh M;
    triplet mass, stiff, dump;
    csr M_csr, K_csr, D_csr;
    SimulationConfig config;
    double *u, *v, *a, *F, *dirichletValues;

    // Control parameters
    double eps1 = 0.00012;
    double eps2 = 0.00004;
    double alpha = 2.0;
    int gravity = FALSE;
    int Damping = TRUE;
    double dt = 0.00006;
    int nSteps = 80000;
    int saveVKT = FALSE;
    int saveData = TRUE;
    // double omega = 165.533616484098; // Old version
    // ---- NEW versions ---
    // double omega = 165.640632154201; // case mesh_6 2280 - beam5
    // double omega = 166.239515265078; // case mesh_5 1108 - beam4
    // double omega = 167.25340512031; // case mesh_4 554 - beam3
    // double omega = 169.453792191287; // case mesh_3 282 - beam2
    // double omega = 173.944119943247; // case mesh_2 141 - beam1
    double omega = 180.254963005736; // case mesh_1 68 - beam

    double f_FFT = 1 / (2 * dt);
    double T = dt * nSteps;

    Mesh_Load(&M, "Elasticity_Geo/Beam2D_bend.msh");
    int nDOF = 2 * M.NPoints;

    printf("FFT frequency: %g\n", f_FFT);
    printf("Size of problem (2*points): %d, %d\n", nDOF, M.NPoints);
    printf("Size of problem (elements): %d \n", M.NElements);
    printf("Signal duration T: %g\n", T);

    InitializeDynamicProblem(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, &u, &v, &a, &F, &dirichletValues, &nDOF, alpha, eps1, eps2, gravity, &config);
    ExecuteDynamicAnalysis(&M, &M_csr, &K_csr, &D_csr, F, u, v, a, nDOF, omega, eps1, eps2, dt, nSteps, Damping, saveVKT, saveData, &config);
    CleanupDynamicProblem(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, u, v, a, F, dirichletValues);

#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef POISSON3D

    mesh M;
    triplet stiff;
    csr A_csr;
    double *b = NULL, *dirichletValues = NULL, *x = NULL;

    // Control parameters
    int solverType = 1;    // 1 = UMFPACK, 2 = Iterative Solver
    double tolerance = 1e-9;

    Mesh_Load(&M, "Elasticity_Geo/3D_test_case_tri_tetra.msh");
    // Mesh_Load(&M, "Elasticity_Geo/3D_test_case_quad_hexa.msh");
    // Mesh_Load(&M, "Elasticity_Geo/3D_test_case_WEDGE_alt.msh");
    
    // Element_vertex_print(&M);
    // Mesh_Classify_Nodes(&M);
    InitializePoissonProblem3D(&M, &stiff, &A_csr, &b, &dirichletValues, &x);
    ExecutePoissonAnalysis3D(&M, &A_csr, b, x, M.NPoints, solverType, tolerance);

    // double L2NormError3D = computeL2Error3D(&M, x, AnalyticalSolution3D);
    // printf("NELEMENTS:%d \t Points:%d\n", M.NElements, M.NPoints);
    // printf("L2 Norm of the solution error: %g\n", L2NormError3D);

    // for (int i = 0; i < M.NPoints; i++)
    // {
    //     printf("b[%d] = %g\n", i, b[i]);
    // }
    
    // printf("\n");

    for (int i = 0; i < M.NPoints; i++)
    {
        printf("x[%d] = %g\n", i, x[i]);
    }

    // for (int i = 0; i < M.NPoints; i++)
    // {
    //     printf("x[%d] (%g, %g, %g) = %g (anlaytic result: %g)\n", i, M.x[i], M.y[i], M.z[i], x[i], sin(M_PI*M.x[i]) * sin(M_PI*M.y[i]) * sin(M_PI*M.z[i]) );
    // }
    // double max_error = 0.0;
    // int index = 0;
    // for (int i = 0; i < M.NPoints; i++)
    // {
        
    //     double error = fabs( x[i] - sin(M_PI*M.x[i]) * sin(M_PI*M.y[i]) * sin(M_PI*M.z[i]) );

    //     printf("Error_node %d (%g, %g, %g) = %g \n", i, M.x[i], M.y[i], M.z[i], error );

    //     if (error > max_error)
    //     {
    //         max_error = error;
    //         index = i;
    //     }
        

    // }

    // printf("Max_error[%d] = %g\n", index, max_error);

    exportNumericalSolutionToVTK_3D(&M, x, "Poisson3D_tetra_DIR_NEU_ROB.vtk" );
    
    CleanupPoissonProblem3D(&stiff, &A_csr, b, dirichletValues, x);

#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef STATICS3D
    mesh M;
    triplet stiff;
    csr K_csr;
    double *b = NULL, *dirichletValues = NULL, *x = NULL;

    // Control parameters
    int gravity = FALSE;
    int solverType = 1;    // 1 = UMFPACK, 2 = Iterative Solver
    double tolerance = 1e-9;

    int nDOF = 3 * M.NPoints;

    double F_magnitude = 100;
    double azimuth_deg = 270.0; //xy - plane
    double elevation_deg = 0.0; // z - direction

    Mesh_Load(&M, "Elasticity_Geo/MESH_3D/BEAM_3D_CANTILEVER.msh");
    // Mesh_Load(&M, "Elasticity_Geo/Vocals.msh");

    Element_vertex_print(&M, FALSE);
    printf("case1\n");
    getchar();
    printf("case2\n");

    InitializeStaticsProblem3D(&M, &stiff, &K_csr, &b, &dirichletValues, &x, F_magnitude, azimuth_deg, elevation_deg, gravity);
    
    ExecuteStaticsAnalysis3D(&M, &K_csr, b, x, nDOF, solverType, tolerance);
    // saveVectorToFile(b, 3*M.NPoints, 1, "Linear_load_vector.txt");
    exportDisplacementToVTK_3D(&M, x, "Sillicone_Rubber_100N_4.vtk");

    double max_disp = 0;
    for (int i = 0; i < 3 * M.NPoints; i++) {
        if (fabs(x[i]) > max_disp) max_disp = fabs(x[i]);
    }
    printf("Max displacement = %.10e\n", max_disp);

    CleanupStaticsProblem3D(&stiff, &K_csr, b, dirichletValues, x);

#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef STATICS3D_ORTHO
mesh M;
    triplet stiff;
    csr K_csr;
    double *b = NULL, *dirichletValues = NULL, *x = NULL;

    // Control parameters
    int gravity = FALSE;
    int solverType = 1;    // 1 = UMFPACK, 2 = Iterative Solver
    double tolerance = 1e-9;

    int nDOF = 3 * M.NPoints;

    double F_magnitude = 1e2;
    double azimuth_deg = 270.0; //xy - plane
    double elevation_deg = 0.0; // z - direction

    double angleX = 0.0;
    double angleY = 0.0;
    double angleZ = 90.0;

    Mesh_Load(&M,"Elasticity_Geo/AGARD_v2.msh");

    InitializeStaticsProblem3D_Orthotropic(&M, &stiff, &K_csr, &b, &dirichletValues, &x, F_magnitude, azimuth_deg, elevation_deg, gravity, angleX, angleY, angleZ);
    
    ExecuteStaticsAnalysis3D_Orthotropic(&M, &K_csr, b, x, nDOF, solverType, tolerance);
    
    // saveVectorToFile(b, 3*M.NPoints, 1, "Linear_load_vector.txt");
    exportDisplacementToVTK_3D(&M, x, "Orthotropic_wing_statics.vtk");

    CleanupStaticsProblem3D_Orthotropic(&stiff, &K_csr, b, dirichletValues, x);
#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef DYNAMICS3D
    mesh M;
    triplet mass, stiff, dump;
    csr M_csr, K_csr, D_csr;
    SimulationConfig config;
    double *u, *v, *a, *F, *dirichletValues;

    // Control parameters
    double eps1 = 0.00012;
    double eps2 = 0.00004;
    double alpha = 2.0;
    int gravity = FALSE;
    int Damping = FALSE;
    double dt = 0.00006;
    int nSteps = 10000;
    int saveVKT = FALSE;
    int saveData = TRUE;
    
    double omega = 2.57129586186301;
    
    // ---- NEW versions ---
    // double omega = 88.5883264487756;//Version cl0.21
    // double omega = 122.02757149261; // Case_1 Transfinite Line = 10
    // double omega = 140.155252310056; // Case_2 Transfinite Line = 8
    // double omega = 173.019395850875; // Case_3 Transfinite Line = 6

    double F_magnitude = 0.0;
    double azimuth_deg = 270.0; //xy - plane
    double elevation_deg = 0.0; // z - direction

    double f_FFT = 1 / (2 * dt);
    double T = dt * nSteps;

    Mesh_Load(&M, "Elasticity_Geo/MESH_3D/BEAM_3D_CANTILEVER.msh");
    int nDOF = 3 * M.NPoints;

    Element_vertex_print(&M, FALSE);
    getchar();
    printf("FFT frequency: %g\n", f_FFT);
    printf("Size of problem (3*points): %d, %d\n", nDOF, M.NPoints);
    printf("Size of problem (elements): %d \n", M.NElements);
    printf("Signal duration T: %g\n", T);
    
    InitializeDynamicProblem3D(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, &u, &v, &a, &F, &dirichletValues, &nDOF, alpha, F_magnitude, azimuth_deg, elevation_deg, eps1, eps2, gravity, &config);
    ExecuteDynamicAnalysis3D(&M, &M_csr, &K_csr, &D_csr, F, u, v, a, nDOF, omega, eps1, eps2, dt, nSteps, Damping, saveVKT, saveData, &config);
    CleanupDynamicProblem3D(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, u, v, a, F, dirichletValues);
#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef DYNAMICS3D_ORTHO
    mesh M;
    triplet mass, stiff, dump;
    csr M_csr, K_csr, D_csr;
    SimulationConfig config;
    double *u, *v, *a, *F, *dirichletValues;

    // Control parameters
    double eps1 = 0.00012;
    double eps2 = 0.000004;
    double alpha = 2.0;
    int gravity = FALSE;
    int Damping = FALSE;
    double dt = 0.00003;
    int nSteps = 80000;
    int saveVKT = TRUE;
    int saveData = TRUE;
    
    // ---- NEW versions ---
    double omega = 94.9994454311559;

    double F_magnitude = 0.0;
    double azimuth_deg = 270.0; //xy - plane
    double elevation_deg = 0.0; // z - direction

    double f_FFT = 1 / (2 * dt);
    double T = dt * nSteps;

    double angleX = 0.0; //YZ
    double angleY = 0.0; //XZ interested in this one
    double angleZ = 0.0; //XY

    Mesh_Load(&M,"Elasticity_Geo/AGARD_v2.msh");
    // Mesh_Load(&M,"Elasticity_Geo/MESH_3D/BEAM_3D_CANTILEVER.msh");
    int nDOF = 3 * M.NPoints;

    Element_vertex_print(&M, FALSE);
    getchar();
    printf("FFT frequency: %g\n", f_FFT);
    printf("Size of problem (3*points): %d, %d\n", nDOF, M.NPoints);
    printf("Size of problem (elements): %d \n", M.NElements);
    printf("Signal duration T: %g\n", T);
    
    InitializeDynamicProblem3D_Orthotropic(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, &u, &v, &a, &F, &dirichletValues, &nDOF, alpha, F_magnitude, azimuth_deg, elevation_deg, eps1, eps2, gravity, angleX, angleY, angleZ, Damping, &config);
    ExecuteDynamicAnalysis3D_Orthotropic(&M, &M_csr, &K_csr, &D_csr, F, u, v, a, nDOF, omega, eps1, eps2, dt, nSteps, Damping, saveVKT, saveData, &config);
    CleanupDynamicProblem3D_Orthotropic(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, Damping, u, v, a, F, dirichletValues);
#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef NONLINEAR
    mesh M;
    triplet stiff;
    csr K_csr;
    double *displacement = NULL, *residual = NULL, *f_ext = NULL, *P_int = NULL, *dirichletValues = NULL;
    
    int maxIter = 2000;
    double tolerance = 1e-16;

    double F_magnitude = 100.0; //WAS 100
    double azimuth_deg = 270.0; //xy - plane
    double elevation_deg = 0.0; // z - direction

    Mesh_Load(&M, "Elasticity_Geo/MESH_3D/BEAM_3D_CANTILEVER.msh");

    Element_vertex_print(&M, FALSE);
    printf("M.NPoint: %d", M.NPoints);
    getchar();
    // printf("case2\n");
    InitializeNonlinearElasticityProblem(&M, &stiff, &K_csr, &displacement, &residual, &f_ext, &P_int, &dirichletValues, F_magnitude, azimuth_deg, elevation_deg);    
    // saveVectorToFile(displacement, 3*M.NPoints, 1, "Nonlinear_load_vector.txt");
    // ExecuteNonlinearElasticityAnalysis3D(&M, stiff, displacement, residual, dirichletValues, F_magnitude, azimuth_deg, elevation_deg, maxIter, tolerance);
    ExecuteNonlinearElasticityAnalysis3D(&M, stiff, displacement, residual, f_ext, P_int, dirichletValues, F_magnitude, azimuth_deg, elevation_deg, maxIter, tolerance);
    // exportDisplacementToVTK_3D(&M, displacement, "NonlinearElasticity_SVK_BUCKLING_SILICONE.vtk");
    exportDisplacementToVTK_3D(&M, displacement, "Sillicone_Rubber_100N_Nonlinear_4.vtk");


    double max_disp = 0;
    for (int i = 0; i < 3 * M.NPoints; i++) {
        if (fabs(displacement[i]) > max_disp) max_disp = fabs(displacement[i]);
    }
    printf("Max displacement = %.10e\n", max_disp);

    CleanupNonlinearElasticityProblem(&M, &stiff, &K_csr, displacement, residual, dirichletValues);
    
#endif
//--------------------------------------------------------------------------------------------------------------
#ifdef NONLINEAR_DYNAMICS
    mesh M;
    triplet mass, stiff, dump;
    csr M_csr, K_csr, D_csr;
    SimulationConfig config;
    double *u = NULL, *v = NULL, *a = NULL, *F = NULL, *dirichletValues = NULL;

    // Control parameters
    double eps1 = 0.00012;
    double eps2 = 0.00004;
    double alpha = 2.0;
    int gravity = FALSE;
    int Damping = FALSE;
    double dt = 0.00006;
    int nSteps = 60;
    int maxNewton = 20;
    int maxIter = 30;
    double tolerance = 1e-8;
    int saveVKT = TRUE;
    int saveData = FALSE;
    
    // ---- NEW versions ---
    // double omega = 88.5883264487756;//Version cl0.21
    // double omega = 122.02757149261; // Case_1 Transfinite Line = 10
    // double omega = 140.155252310056; // Case_2 Transfinite Line = 8
    // double omega = 173.019395850875; // Case_3 Transfinite Line = 6

    double F_magnitude = 100;
    double azimuth_deg = 270.0; //xy - plane
    double elevation_deg = 0.0; // z - direction

    double f_FFT = 1 / (2 * dt);
    double T = dt * nSteps;

    Mesh_Load(&M, "Elasticity_Geo/MESH_3D/BEAM_3D_CANTILEVER.msh");
    int nDOF = 3 * M.NPoints;

    Element_vertex_print(&M, FALSE);
    getchar();
    printf("FFT frequency: %g\n", f_FFT);
    printf("Size of problem (3*points): %d, %d\n", nDOF, M.NPoints);
    printf("Size of problem (elements): %d \n", M.NElements);
    printf("Signal duration T: %g\n", T);

    InitializeNonlinearElasticityProblem_Dynamics(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, &u, &v, &a, &F, &dirichletValues, maxIter, tolerance, maxNewton, dt, nSteps, Damping, saveVKT, saveData, &nDOF, alpha, F_magnitude, azimuth_deg, elevation_deg, eps1, eps2, gravity, &config);
    
    // InitializeNonlinearElasticityProblem_Dynamics(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, &u, &v, &a, &F, &dirichletValues, 20, 1e-8, 20, dt, nSteps, Damping, saveVKT, saveData, &nDOF, alpha, F_magnitude, azimuth_deg, elevation_deg, eps1, eps2, gravity, &config);
    CleanupNonlinearElasticityProblem_Dynamics(&M, &mass, &stiff, &dump, &M_csr, &K_csr, &D_csr, u, v, a, F, dirichletValues);
#endif
    return 0; 
}

//--------------------------------------------------------------------------------------------------------------

void Element_vertex_print(mesh *M, int Vertices)
{
    printf("\n");
    printf("NElemeent: %d\n", M->NElements);
    printf("NTriangles: %d\n", M->NTriangles);
    printf("NQuads: %d\n", M->NQuads);
    printf("NTetrahedra: %d\n", M->NTetrahedra);
    printf("NHexahedra: %d\n", M->NHexahedra);
    printf("NPyramid: %d\n", M->NPyramids);
    printf("NWedge: %d\n", M->NWedges);
    printf("\n");

    if (Vertices == TRUE)
    {    
        printf("\nMesh Elements:\n");

        // Offsets for different element types
        int startQuad = 3 * M->NTriangles;                            // Quadrangles start
        int startTetra = startQuad + 4 * M->NQuads;                   // Tetrahedra start
        int startHex = startTetra + 4 * M->NTetrahedra;               // Hexahedra start
        int startPrism = startHex + 8 * M->NHexahedra;                // Prisms start
        int startPyramid = startPrism + 5 * M->NWedges;               // Pyramids start

        for (int i = 0; i < M->NElements; i++) {
            if (i < M->NTriangles) {
                // Print Triangle (3 nodes)
                printf("Triangle %d: ", i + 1);
                for (int j = 0; j < 3; j++) {
                    int vertex = M->ElementListOfVertices[3 * i + j];
                    printf("%d ", vertex + 1);  // 1-based indexing
                }
                printf("\n");
            } 
            else if (i >= M->NTriangles && i < (M->NTriangles + M->NQuads)) {
                // Print Quadrangle (4 nodes)
                int quadIndex = i - M->NTriangles;
                printf("Quadrangle %d: ", quadIndex + 1);
                for (int j = 0; j < 4; j++) {
                    int vertex = M->ElementListOfVertices[startQuad + 4 * quadIndex + j];
                    printf("%d ", vertex + 1);  // 1-based indexing
                }
                printf("\n");
            } 
            else if (i >= (M->NTriangles + M->NQuads) && i < (M->NTriangles + M->NQuads + M->NTetrahedra)) {
                // Print Tetrahedron (4 nodes)
                int tetraIndex = i - (M->NTriangles + M->NQuads);
                printf("Tetrahedron %d: ", tetraIndex + 1);
                for (int j = 0; j < 4; j++) {
                    int vertex = M->ElementListOfVertices[startTetra + 4 * tetraIndex + j];
                    printf("%d ", vertex + 1);  // 1-based indexing
                }
                printf("\n");
            } 
            else if (i >= (M->NTriangles + M->NQuads + M->NTetrahedra) && i < (M->NTriangles + M->NQuads + M->NTetrahedra + M->NHexahedra)) {
                // Print Hexahedron (8 nodes)
                int hexIndex = i - (M->NTriangles + M->NQuads + M->NTetrahedra);
                printf("Hexahedron %d: ", hexIndex + 1);
                for (int j = 0; j < 8; j++) {
                    int vertex = M->ElementListOfVertices[startHex + 8 * hexIndex + j];
                    printf("%d ", vertex + 1);  // 1-based indexing
                }
                printf("\n");
            } 
            else if (i >= (M->NTriangles + M->NQuads + M->NTetrahedra + M->NHexahedra) && i < (M->NTriangles + M->NQuads + M->NTetrahedra + M->NHexahedra + M->NWedges)) {
                // Print Prism (6 nodes)
                int prismIndex = i - (M->NTriangles + M->NQuads + M->NTetrahedra + M->NHexahedra);
                printf("Prism (Wedge) %d: ", prismIndex + 1);
                for (int j = 0; j < 6; j++) {
                    int vertex = M->ElementListOfVertices[startPrism + 6 * prismIndex + j];
                    printf("%d ", vertex + 1);  // 1-based indexing
                }
                printf("\n");
            } 
            else if (i >= (M->NTriangles + M->NQuads + M->NTetrahedra + M->NHexahedra + M->NWedges) && i < M->NElements) {
                // Print Pyramid (5 nodes)
                int pyramidIndex = i - (M->NTriangles + M->NQuads + M->NTetrahedra + M->NHexahedra + M->NWedges);
                printf("Pyramid %d: ", pyramidIndex + 1);
                for (int j = 0; j < 5; j++) {
                    int vertex = M->ElementListOfVertices[startPyramid + 5 * pyramidIndex + j];
                    printf("%d ", vertex + 1);  // 1-based indexing
                }
                printf("\n");
            } 
            else {
                printf("Unknown element at index %d\n", i + 1);
            }
        }
        printf("\n");
    }
}

void Mesh_Node_Coordinates_Print(mesh *M)
{
    printf("\nMesh Node Coordinates:\n");

    for (int i = 0; i < M->NPoints; i++) {
        printf("Node %d: x = %f, y = %f",
               i + 1, // 1-based indexing for node ID
               M->x[i], // x-coordinate
               M->y[i]); // y-coordinate

        if (M->z != NULL) { // Check if the mesh has z-coordinates for 3D
            printf(", z = %f", M->z[i]); // z-coordinate
        }

        printf("\n");
    }
    printf("\n");
}

void Mesh_Classify_Nodes(mesh *M)
{
    // Allocate memory to track surface nodes
    bool *isSurfaceNode = calloc(M->NPoints, sizeof(bool));

    // Traverse through all surface elements (Triangles and Quadrangles)
    for (int i = 0; i < M->NTriangles; i++) {
        isSurfaceNode[M->ElementListOfVertices[3 * i + 0]] = true;
        isSurfaceNode[M->ElementListOfVertices[3 * i + 1]] = true;
        isSurfaceNode[M->ElementListOfVertices[3 * i + 2]] = true;
    }

    int startQuad = 3 * M->NTriangles;
    for (int i = 0; i < M->NQuads; i++) {
        isSurfaceNode[M->ElementListOfVertices[startQuad + 4 * i + 0]] = true;
        isSurfaceNode[M->ElementListOfVertices[startQuad + 4 * i + 1]] = true;
        isSurfaceNode[M->ElementListOfVertices[startQuad + 4 * i + 2]] = true;
        isSurfaceNode[M->ElementListOfVertices[startQuad + 4 * i + 3]] = true;
    }

    // Output which nodes are on the surface and which are internal
    printf("\nNode Classification:\n");
    for (int i = 0; i < M->NPoints; i++) {
        if (isSurfaceNode[i]) {
            printf("Node %d (%.3f, %.3f, %.3f) is a SURFACE node.\n", 
                   i + 1, M->x[i], M->y[i], M->z[i]);
        } else {
            printf("Node %d (%.3f, %.3f, %.3f) is an INTERNAL node.\n", 
                   i + 1, M->x[i], M->y[i], M->z[i]);
        }
    }

    // Cleanup allocated memory
    free(isSurfaceNode);
}

