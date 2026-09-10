#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "problems.h"
#include "mesh.h"
#include "gelement3D.h"
#include "Nonlinear.h"
#include "saving.h"
#include "constants.h"
/* ------------------------------------------------------------------------------------------------- */

void InitializePoissonProblem(mesh *M, triplet *stiff, csr *A_csr, double **b, double **dirichletValues, double **x, IntegrationType integrationType) {
    int nDOF = M->NPoints;

    // Allocate triplet, load vector, and solution vector
    Triplet_Allocate(stiff, nDOF, M->NElements * 100);
    *b = calloc(nDOF, sizeof(double));
    *dirichletValues = calloc(nDOF, sizeof(double));
    *x = calloc(nDOF, sizeof(double)); // Allocate and initialize solution vector

    // Set Dirichlet boundary conditions
    for (int i = 0; i < M->NPoints; i++) {
        if (M->PointMark[i] == DIRICHLET_BC) {
            (*dirichletValues)[i] = sin(M_PI * M->x[i]) * sin(M_PI * M->y[i]);
        }
    }

    Mesh_Discretize_Poisson(M, f, stiff, *b, *dirichletValues, psi, ur, integrationType);

    CSR_InitFromTriplet(A_csr, stiff);

    saveSparseMatrixToFile(A_csr, "Poisson_Stiffness_Matrix.dat");
}

void ExecutePoissonAnalysis(mesh *M, csr *A_csr, double *b, double *x, int nDOF, int solverType, double tolerance) {
    printf("Executing Poisson analysis...\n");

    if (solverType == 1) {
        CSR_Solve(A_csr, b, x);
    } else {
        ConjugateGradient(A_csr, b, x, 1e5, tolerance);
    }

    saveVectorToFile(x, nDOF, 1, "Poisson_Displacement_Vector.dat");

    printf("Poisson analysis completed.\n");
}

void CleanupPoissonProblem(triplet *stiff, csr *A_csr, double *b, double *dirichletValues, double *x) {
    Triplet_Free(stiff);
    CSR_Free(A_csr);
    free(b);
    free(dirichletValues);
    free(x);
    printf("Poisson problem resources cleaned up.\n");
}

/* ------------------------------------------------------------------------------------------------- */

void InitializePoissonProblem3D(mesh *M, triplet *stiff, csr *A_csr, double **b, double **dirichletValues, double **x) {
    int nDOF = M->NPoints;

    Triplet_Allocate(stiff, nDOF, M->NElements * 1000);
    *b = calloc(nDOF, sizeof(double));
    *dirichletValues = calloc(nDOF, sizeof(double));
    *x = calloc(nDOF, sizeof(double));

    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            // printf("Element %d: %d\n", elem, K.type);
            for (int i = 0; i < K.Nvertex; i++)
            {
                int node = K.idxNode[i];
                (*dirichletValues)[node] = 1.0 * ( sin(M_PI * M->x[node]) * sin(M_PI * M->y[node]) * sin(M_PI * M->z[node]) );
                // printf("Node: %d\n", node);
            }
            // printf("\n");
        }
    }

    Mesh_Discretize_Poisson_3D(M, f3D, stiff, *b, *dirichletValues, psi3D, ur3D);

    CSR_InitFromTriplet(A_csr, stiff);

    saveSparseMatrixToFile(A_csr, "Poisson_3D_Stiffness_Matrix.dat");
}

void ExecutePoissonAnalysis3D(mesh *M, csr *A_csr, double *b, double *x, int nDOF, int solverType, double tolerance) {
    printf("Executing 3D Poisson analysis...\n");

    if (solverType == 1) {
        CSR_Solve(A_csr, b, x);
    } else {
        ConjugateGradient(A_csr, b, x, 1e5, tolerance);
    }

    saveVectorToFile(x, nDOF, 1, "Poisson_3D_Displacement_Vector.dat");

    printf("3D Poisson analysis completed.\n");
}

void CleanupPoissonProblem3D(triplet *stiff, csr *A_csr, double *b, double *dirichletValues, double *x) {
    Triplet_Free(stiff);
    CSR_Free(A_csr);
    free(b);
    free(dirichletValues);
    free(x);
    printf("3D Poisson problem resources cleaned up.\n");
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeStaticsProblem(mesh *M, triplet *stiff, csr *K_csr, double **b, double **dirichletValues, double **x, double angle_degrees, int gravity) {
    int nDOF = 2 * M->NPoints;

    // Allocate stiffness triplet, load vector, and solution vector
    Triplet_Allocate(stiff, nDOF, M->NElements * 100);
    *b = calloc(nDOF, sizeof(double));
    *dirichletValues = calloc(nDOF, sizeof(double));
    *x = calloc(nDOF, sizeof(double)); // Allocate and initialize solution vector

    // Set Dirichlet boundary conditions
    for (int i = 0; i < M->NPoints; i++) {
        if (M->PointMark[i] == DIRICHLET_BC) {
            (*dirichletValues)[2 * i] = 0.0;
            (*dirichletValues)[2 * i + 1] = 0.0;
        }
    }

    // Discretize the statics problem
    Mesh_Discretize_Statics(M, *b, stiff, *dirichletValues, psi_force, angle_degrees, gravity);

    // Convert triplet to CSR format
    CSR_InitFromTriplet(K_csr, stiff);

    // Save stiffness matrix
    // saveSparseMatrixToFile(K_csr, "Statics_Stiffness_Matrix.dat");
}

void ExecuteStaticsAnalysis(mesh *M, csr *K_csr, double *b, double *x, int nDOF, int solverType, double tolerance) {
    printf("Executing statics analysis...\n");

    // Solve the system K * x = b
    if (solverType == 1) {
        // Direct solver (UMFPACK)
        CSR_Solve(K_csr, b, x);
    } else {
        // Iterative solver (e.g., Conjugate Gradient)
        ConjugateGradient(K_csr, b, x, 1e5, tolerance);
    }

    // Save the solution
    saveVectorToFile(x, nDOF, 1, "Statics_Displacement_Vector.dat");

    printf("Statics analysis completed.\n");
}

void CleanupStaticsProblem(triplet *stiff, csr *K_csr, double *b, double *dirichletValues, double *x) {
    Triplet_Free(stiff);
    CSR_Free(K_csr);
    free(b);
    free(dirichletValues);
    free(x);
    printf("Statics problem resources cleaned up.\n");
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeDynamicProblem(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, int *nDOF, double alpha, double eps1, double eps2, int gravity, SimulationConfig *config) {
    
    initializeSimulationConfig(config);

    Triplet_Allocate(stiff, *nDOF, M->NElements * 100);
    Triplet_Allocate(mass, *nDOF, M->NElements * 100);
    Triplet_Allocate(dump, *nDOF, M->NElements * 100);

    *u = calloc(*nDOF, sizeof(double));
    *v = calloc(*nDOF, sizeof(double));
    *a = calloc(*nDOF, sizeof(double));
    *F = calloc(*nDOF, sizeof(double));
    *dirichletValues = calloc(*nDOF, sizeof(double));

    for (int i = 0; i < M->NPoints; i++) {
        if (M->PointMark[i] == DIRICHLET_BC) {
            (*dirichletValues)[2 * i] = 0.0;
            (*dirichletValues)[2 * i + 1] = 0.0;
        }
    }

    for (int i = 0; i < M->NPoints; i++) {
        int dof_index = 2 * i;
        double x = M->x[i];
        (*v)[dof_index] = 0.0;
        (*v)[dof_index + 1] = - 2*x;
    }

    // for (int i = 0; i < M->NPoints; i++) {
    //     int dof_index = 2 * i;
    //     // double x = M->x[i];
    //     (*F)[dof_index] = 0.0;
    //     (*F)[dof_index + 1] = 0.0;//- 2*x;
    // }
    
    Mesh_Discretize_Dynamics(M, *F, mass, stiff, dump, *dirichletValues, psi_force, 180.0, eps1, eps2, gravity);

    CSR_InitFromTriplet(M_csr, mass);
    CSR_InitFromTriplet(K_csr, stiff);
    CSR_InitFromTriplet(D_csr, dump);

    saveSparseMatrixToFile(M_csr, "Matrices_modal/NUnscaled_Dyn_Mass_damp_testing68.dat"); //_damp
    saveSparseMatrixToFile(K_csr, "Matrices_modal/NUnscaled_Dyn_Stiff_damp_testing68.dat");
    saveSparseMatrixToFile(D_csr, "Matrices_modal/NUnscaled_Dyn_Dump_damp_testing68.dat");
    printf("Stop after unscaled matrices loaded:\n");
    getchar();
}

void ExecuteDynamicAnalysis(mesh *M, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, int nDOF, double omega, double eps1, double eps2, double dt, int nSteps, int Damping, int saveVKT, int saveData, SimulationConfig *config) {
    
    scaleMatricesForDirichlet(M, M_csr, Damping ? D_csr : NULL, omega, eps1, eps2, 2, Damping);

    saveSparseMatrixToFile(M_csr, "Matrices_modal/NEW_Scaled_Dyn_Mass_damp_testing68.dat");
    saveSparseMatrixToFile(K_csr, "Matrices_modal/NEW_Scaled_Dyn_Stiff_damp_testing68.dat");
    saveSparseMatrixToFile(D_csr, "Matrices_modal/NEW_Scaled_Dyn_Dumpda_damp_testing68.dat");

    NewmarkIntegrate(M, M_csr, K_csr, D_csr, F, u, v, a, dt, nSteps, Damping, saveVKT, saveData, config);
}

void CleanupDynamicProblem(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double *u, double *v, double *a, double *F, double *dirichletValues) {
    Triplet_Free(mass);
    Triplet_Free(stiff);
    Triplet_Free(dump);

    CSR_Free(M_csr);
    CSR_Free(K_csr);
    CSR_Free(D_csr);

    free(u);
    free(v);
    free(a);
    free(F);
    free(dirichletValues);

    Mesh_Free(M);
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeStaticsProblem3D(mesh *M, triplet *stiff, csr *K_csr, double **b, double **dirichletValues, double **x, double F_magnitude, double azimuth_deg, double elevation_deg, int gravity) {
    int nDOF = 3 * M->NPoints;

    int non_zero_allocated = (int)( 1.4 * ( M->NTetrahedra * 144 + M->NHexahedra * 576 + M->NWedges * 324 + M->NPyramids * 225 ) );
    Triplet_Allocate(stiff, nDOF, non_zero_allocated);
    *b = calloc(nDOF, sizeof(double));
    *dirichletValues = calloc(nDOF, sizeof(double));
    *x = calloc(nDOF, sizeof(double)); // Allocate and initialize solution vector

    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            for (int i = 0; i < K.Nvertex; i++) {
                int node = K.idxNode[i];
                // Set ALL vector components for this node
                (*dirichletValues)[3 * node    ] = 0.0;
                (*dirichletValues)[3 * node + 1] = 0.0;
                (*dirichletValues)[3 * node + 2] = 0.0;
            }
        }
    }
    
    double angleX = 0.0;
    double angleY = 0.0;
    double angleZ = 0.0;
    
    Mesh_Discretize_Statics3D(M, *b, stiff, *dirichletValues, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg, gravity, angleX, angleY, angleZ);

    // printf("Mesh_discretize_3D (F1)\n");
    // Convert triplet to CSR format
    CSR_InitFromTriplet(K_csr, stiff);
    // printf("CSR_from_triplet (F1)\n");
    // Save stiffness matrix
    // saveSparseMatrixToFile(K_csr, "Statics_Stiffness_Matrix.dat");
}

void ExecuteStaticsAnalysis3D(mesh *M, csr *K_csr, double *b, double *x, int nDOF, int solverType, double tolerance) {
    printf("Executing statics analysis...\n");

    if (solverType == 1) {
        CSR_Solve(K_csr, b, x);
    } else {
        
        ConjugateGradient(K_csr, b, x, 1e5, tolerance);
    }

    // saveVectorToFile(x, nDOF, 1, "Statics_Displacement_Vector.dat");

    printf("Statics analysis completed.\n");
}

void CleanupStaticsProblem3D(triplet *stiff, csr *K_csr, double *b, double *dirichletValues, double *x) {
    Triplet_Free(stiff);
    CSR_Free(K_csr);
    free(b);
    free(dirichletValues);
    free(x);
    printf("Statics 3D problem resources cleaned up.\n");
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeStaticsProblem3D_Orthotropic(mesh *M, triplet *stiff, csr *K_csr, double **b, double **dirichletValues, double **x, double F_magnitude, double azimuth_deg, double elevation_deg, int gravity, double angleX, double angleY, double angleZ) {
    int nDOF = 3 * M->NPoints;

    // Allocate memory for triplet and vectors
    int non_zero_allocated = (int)(600.4 * (M->NTetrahedra * 144 + M->NHexahedra * 576 + M->NWedges * 324 + M->NPyramids * 225));
    Triplet_Allocate(stiff, nDOF, non_zero_allocated);
    *b = calloc(nDOF, sizeof(double));
    *dirichletValues = calloc(nDOF, sizeof(double));
    *x = calloc(nDOF, sizeof(double)); // Allocate and initialize solution vector
    printf("Variables initialized (Orthotropic)\n");

    // Set Dirichlet boundary conditions
    // for (int i = 0; i < M->NPoints; i++) {
    //     if (M->PointMark[i] == DIRICHLET_BC) {
    //         (*dirichletValues)[3 * i] = 0.0;
    //         (*dirichletValues)[3 * i + 1] = 0.0;
    //         (*dirichletValues)[3 * i + 2] = 0.0;
    //     }
    // }
    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            for (int i = 0; i < K.Nvertex; i++) {
                int node = K.idxNode[i];
                // Set ALL vector components for this node
                (*dirichletValues)[3 * node    ] = 0.0;
                (*dirichletValues)[3 * node + 1] = 0.0;
                (*dirichletValues)[3 * node + 2] = 0.0;
            }
        }
    }
    printf("Dirichlet initialized (Orthotropic)\n");

    // Discretize the statics problem with orthotropic materials
    Mesh_Discretize_Statics3D_Orthotropic(M, *b, stiff, *dirichletValues, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg, gravity, angleX, angleY, angleZ);
    printf("Mesh_discretize_3D_Orthotropic (Orthotropic)\n");

    // Convert triplet to CSR format
    CSR_InitFromTriplet(K_csr, stiff);
    printf("CSR_from_triplet (Orthotropic)\n");

    // saveSparseMatrixToFile(K_csr, "Orthotropic_Stiffness_Matrix.dat");
}

void ExecuteStaticsAnalysis3D_Orthotropic(mesh *M, csr *K_csr, double *b, double *x, int nDOF, int solverType, double tolerance) {
    printf("Executing statics analysis for orthotropic materials...\n");

    if (solverType == 1) {
        CSR_Solve(K_csr, b, x);
    } else {
        ConjugateGradient(K_csr, b, x, 1e5, tolerance);
    }

    // Optionally, save displacement vector to a file
    // saveVectorToFile(x, nDOF, 1, "Orthotropic_Statics_Displacement_Vector.dat");

    printf("Statics analysis completed for orthotropic materials.\n");
}

void CleanupStaticsProblem3D_Orthotropic(triplet *stiff, csr *K_csr, double *b, double *dirichletValues, double *x) {
    // Free allocated memory
    Triplet_Free(stiff);
    CSR_Free(K_csr);
    free(b);
    free(dirichletValues);
    free(x);
    printf("Orthotropic 3D statics problem resources cleaned up.\n");
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeDynamicProblem3D(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, int *nDOF, double alpha, double F_magnitude, double azimuth_deg, double elevation_deg, double eps1, double eps2, int gravity, SimulationConfig *config) {    

    initializeSimulationConfig(config);

    int non_zero_allocated = (int)( 4.4 * ( M->NTetrahedra * 144 + M->NHexahedra * 576 + M->NWedges * 324 + M->NPyramids * 225 ) );

    Triplet_Allocate(stiff, *nDOF, non_zero_allocated);
    Triplet_Allocate(mass, *nDOF, non_zero_allocated);
    Triplet_Allocate(dump, *nDOF, non_zero_allocated);

    *u = calloc(*nDOF, sizeof(double));
    *v = calloc(*nDOF, sizeof(double));
    *a = calloc(*nDOF, sizeof(double));
    *F = calloc(*nDOF, sizeof(double));
    *dirichletValues = calloc(*nDOF, sizeof(double));

    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            for (int i = 0; i < K.Nvertex; i++) {
                int node = K.idxNode[i];
                // Set ALL vector components for this node
                (*dirichletValues)[3 * node    ] = 0.0;
                (*dirichletValues)[3 * node + 1] = 0.0;
                (*dirichletValues)[3 * node + 2] = 0.0;
            }
        }
    }

    for (int i = 0; i < M->NPoints; i++) {
        int dof_index = 3 * i;
        double x = M->x[i];
        (*v)[dof_index] = 0.0;
        (*v)[dof_index + 1] = -0.5;
        (*v)[dof_index + 2] = 0.0;
    }

    // for (int i = 0; i < M->NPoints; i++) {
    //     int dof_index = 3 * i;
    //     // double x = M->x[i];
    //     (*F)[dof_index] = 0.0;
    //     (*F)[dof_index + 1] = 0.0;//- 2*x;
    //     (*F)[dof_index + 2] = 0.0;//- 2*x;
    // }

    double angleX = 0.0;
    double angleY = 0.0;
    double angleZ = 0.0;

    Mesh_Discretize_Dynamics3D(M, *F, mass, stiff, dump, *dirichletValues, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg, eps1, eps2, gravity, angleX, angleY, angleZ);
    printf("Mesh_Discretize_Dynamics3D done\n");

    CSR_InitFromTriplet(M_csr, mass); printf("CSR_InitFromTriplet MASS done\n");
    CSR_InitFromTriplet(K_csr, stiff); printf("CSR_InitFromTriplet STIFF done\n");
    CSR_InitFromTriplet(D_csr, dump); printf("CSR_InitFromTriplet DAMP done\n");

    saveSparseMatrixToFile(M_csr, "Matrices_modal/Mass_matrix_damping_100K.dat");
    saveSparseMatrixToFile(K_csr, "Matrices_modal/Stiff_matrix_damping_100K.dat");
    saveSparseMatrixToFile(D_csr, "Matrices_modal/Damping_matrix_damping_100K.dat");
    printf("Stop after unscaled matrices loaded:\n");
    getchar();
}

void ExecuteDynamicAnalysis3D(mesh *M, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, int nDOF, double omega, double eps1, double eps2, double dt, int nSteps, int Damping, int saveVKT, int saveData, SimulationConfig *config) {
    
    scaleMatricesForDirichlet(M, M_csr, Damping ? D_csr : NULL, omega, eps1, eps2, 3, Damping);

    saveSparseMatrixToFile(M_csr, "Matrices_modal/Mass_matrix_scaled_damping_100K.dat");
    saveSparseMatrixToFile(K_csr, "Matrices_modal/Stiff_matrix_scaled_damping_100K.dat");
    saveSparseMatrixToFile(D_csr, "Matrices_modal/Damping_matrix_scaled_damping_100K.dat");

    NewmarkIntegrate(M, M_csr, K_csr, D_csr, F, u, v, a, dt, nSteps, Damping, saveVKT, saveData, config);
}

void CleanupDynamicProblem3D(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double *u, double *v, double *a, double *F, double *dirichletValues) {
    Triplet_Free(mass);
    Triplet_Free(stiff);
    Triplet_Free(dump);

    CSR_Free(M_csr);
    CSR_Free(K_csr);
    CSR_Free(D_csr);

    free(u);
    free(v);
    free(a);
    free(F);
    free(dirichletValues);

    Mesh_Free(M);
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeDynamicProblem3D_Orthotropic(mesh *M, triplet *mass, triplet *stiff, triplet *damp, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, int *nDOF, double alpha, double F_magnitude, double azimuth_deg, double elevation_deg, double eps1, double eps2, int gravity, double angleX, double angleY, double angleZ, int Damping, SimulationConfig *config) {    

    initializeSimulationConfig(config);

    int non_zero_allocated = (int)( 4.4 * ( M->NTetrahedra * 144 + M->NHexahedra * 576 + M->NWedges * 324 + M->NPyramids * 225 ) );

    Triplet_Allocate(stiff, *nDOF, non_zero_allocated);
    Triplet_Allocate(mass, *nDOF, non_zero_allocated);
    Triplet_Allocate(damp, *nDOF, non_zero_allocated);

    *u = calloc(*nDOF, sizeof(double));
    *v = calloc(*nDOF, sizeof(double));
    *a = calloc(*nDOF, sizeof(double));
    *F = calloc(*nDOF, sizeof(double));
    *dirichletValues = calloc(*nDOF, sizeof(double));

    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            for (int i = 0; i < K.Nvertex; i++) {
                int node = K.idxNode[i];
                // Set ALL vector components for this node
                (*dirichletValues)[3 * node    ] = 0.0;
                (*dirichletValues)[3 * node + 1] = 0.0;
                (*dirichletValues)[3 * node + 2] = 0.0;
            }
        }
    }

    for (int elem = M->NTriangles + M->NQuads; elem < M->NElements; elem++) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            for (int i = 0; i < K.Nvertex; i++) {
                int node = K.idxNode[i];
                int dof_index = 3 * node; // Use global node number!
                // double y = M->y[node];    // Use global node number!
                (*v)[dof_index]     = 0.0;
                (*v)[dof_index + 1] = -2.0;
                (*v)[dof_index + 2] = 0.0;
            }
    }

    // for (int i = 0; i < M->NPoints; i++) {
    //     int dof_index = 3 * i;
    //     // double x = M->x[i];
    //     (*F)[dof_index] = 0.0;
    //     (*F)[dof_index + 1] = 0.0;//- 2*x;
    //     (*F)[dof_index + 2] = 0.0;//- 2*x;
    // }

    Mesh_Discretize_Dynamics3D_Orthotropic(M, *F, mass, stiff, damp, *dirichletValues, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg, eps1, eps2, gravity, angleX, angleY, angleZ);
    printf("Mesh_Discretize_Dynamics3D done\n");

    CSR_InitFromTriplet(M_csr, mass); printf("CSR_InitFromTriplet MASS done\n");
    CSR_InitFromTriplet(K_csr, stiff); printf("CSR_InitFromTriplet STIFF done\n");

    if (Damping) {
        CSR_InitFromTriplet(D_csr, damp); printf("CSR_InitFromTriplet DAMP done\n");
    }

    // saveSparseMatrixToFile(M_csr, "Wing_matrix_etc/Orthotropic_xyz_0(-90)0_dynamics.dat");
    // saveSparseMatrixToFile(K_csr, "Wing_matrix_etc/Orthotropic_xyz_0(-90)0_stiffness.dat");
    // saveSparseMatrixToFile(D_csr, "Wing_matrix_etc/Orthotropic_xyz_0(-90)0_damping.dat");

    saveSparseMatrixToFile(M_csr, "Wing_testing/Ortotropic_wing_0_mass.dat");
    saveSparseMatrixToFile(K_csr, "Wing_testing/Ortotropic_wing_0_stif.dat");
    if (Damping) {
        saveSparseMatrixToFile(D_csr, "Wing_testing/Ortotropic_wing_0_damping.dat");
    }
    printf("Stop after unscaled matrices loaded:\n");
    getchar();
}

void ExecuteDynamicAnalysis3D_Orthotropic(mesh *M, csr *M_csr, csr *K_csr, csr *D_csr, double *F, double *u, double *v, double *a, int nDOF, double omega, double eps1, double eps2, double dt, int nSteps, int Damping, int saveVKT, int saveData, SimulationConfig *config) {
    
    scaleMatricesForDirichlet(M, M_csr, Damping ? D_csr : NULL, omega, eps1, eps2, 3, Damping);

    // saveSparseMatrixToFile(M_csr, "Wing_matrix_etc/Orthotropic_xyz_0(-90)0_dynamics_scaled.dat");
    // saveSparseMatrixToFile(K_csr, "Wing_matrix_etc/Orthotropic_xyz_0(-90)0_stiffness_scaled.dat");
    // saveSparseMatrixToFile(D_csr, "Wing_matrix_etc/Orthotropic_xyz_0(-90)0_damping_scaled.dat");

    saveSparseMatrixToFile(M_csr, "Wing_testing/Ortotropic_wing_0_mass_scaled.dat");
    saveSparseMatrixToFile(K_csr, "Wing_testing/Ortotropic_wing_0_stiff_scaled.dat");
    if (Damping) {
        saveSparseMatrixToFile(D_csr, "Wing_testing/Ortotropic_wing_0_damping_scaled.dat");
    }
    NewmarkIntegrate(M, M_csr, K_csr, D_csr, F, u, v, a, dt, nSteps, Damping, saveVKT, saveData, config);
}

void CleanupDynamicProblem3D_Orthotropic(mesh *M, triplet *mass, triplet *stiff, triplet *damp, csr *M_csr, csr *K_csr, csr *D_csr, int Damping, double *u, double *v, double *a, double *F, double *dirichletValues) {
    Triplet_Free(mass);
    Triplet_Free(stiff);
    Triplet_Free(damp);

    CSR_Free(M_csr);
    CSR_Free(K_csr);
    if (Damping) {
        CSR_Free(D_csr);
    }

    free(u);
    free(v);
    free(a);
    free(F);
    free(dirichletValues);

    Mesh_Free(M);
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeNonlinearElasticityProblem(mesh *M, triplet *stiff, csr *K_csr, double **displacement, double **residual, double **f_ext, double **P_int, double **dirichletValues, double F_magnitude, double azimuth_deg, double elevation_deg) {
    int nDOF = M->NPoints * M->dim;

    // Allocate necessary triplets and vectors
    int non_zero_allocated = (int)( 2.4 * ( M->NTetrahedra * 576 + M->NHexahedra * 576 + M->NWedges * 576 + M->NPyramids * 576 ) );
    Triplet_Allocate(stiff, nDOF, non_zero_allocated);
    *displacement = calloc(nDOF, sizeof(double));
    *residual = calloc(nDOF, sizeof(double));
    *dirichletValues = calloc(nDOF, sizeof(double));
    *f_ext = calloc(nDOF, sizeof(double));
    *P_int = calloc(nDOF, sizeof(double));

    // for (int i = 0; i < nDOF; i++) {
    //     (*displacement)[i] = -0.5;
    // }
    

    printf("Allocation done\n");

    // Set Dirichlet boundary conditions (example: fixing the first and last points)
    // for (int i = 0; i < M->NPoints; i++) {
    //     if (M->PointMark[i] == DIRICHLET_BC) {
    //         (*dirichletValues)[3 * i] = 0.0;  // Set the displacement for Dirichlet BCs
    //         (*dirichletValues)[3 * i + 1] = 0.0;  // Set the displacement for Dirichlet BCs
    //         (*dirichletValues)[3 * i + 2] = 0.0;  // Set the displacement for Dirichlet BCs
    //     }
    // }
    for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
        if (M->ElementMark[elem] == DIRICHLET_BC) {
            gelement3D K;
            GetElement3D(M, &K, elem);
            for (int i = 0; i < K.Nvertex; i++) {
                int node = K.idxNode[i];
                (*dirichletValues)[3 * node    ] = 0.0;
                (*dirichletValues)[3 * node + 1] = 0.0;
                (*dirichletValues)[3 * node + 2] = 0.0;
            }
        }
    }
}

void ExecuteNonlinearElasticityAnalysis3D(mesh *M, triplet stiff, double *displacement, double *residual, double *f_ext, double *P_int, double *dirichletValues, double F_magnitude, double azimuth_deg, double elevation_deg, int maxIter, double tolerance) {

    printf("Executing nonlinear elasticity analysis...\n");

    // runNonlinearWithLoadSteps(M, displacement,&stiff, residual, f_ext, P_int, dirichletValues, maxIter, tolerance, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg, 100);
    solveNonLinear(M, displacement, &stiff, residual, f_ext, P_int, dirichletValues, maxIter, tolerance, psi_force_3D, F_magnitude, azimuth_deg, elevation_deg);

    printf("Nonlinear elasticity analysis completed.\n");
}

void CleanupNonlinearElasticityProblem(mesh *M, triplet *stiff, csr *K_csr, double *displacement, double *residual, double *dirichletValues) {
    Triplet_Free(stiff);
    free(displacement);
    free(residual);
    free(dirichletValues);
    Mesh_Free(M);
    printf("Nonlinear Elasticity problem resources cleaned up.\n");
}

/* ------------------------------------------------------------------------------------------------- */

void InitializeNonlinearElasticityProblem_Dynamics(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double **u, double **v, double **a, double **F, double **dirichletValues, double maxIter, double tolerance, double maxNewton, double dt, int nSteps, int useDamping, int saveVibration, int saveData, int *nDOF, double alpha, double F_magnitude, double azimuth_deg, double elevation_deg, double eps1, double eps2, int gravity, SimulationConfig *config) {
    initializeSimulationConfig(config);

    int non_zero_allocated = (int)( 4.4 * ( M->NTetrahedra * 144 + M->NHexahedra * 576 + M->NWedges * 324 + M->NPyramids * 225 ) );

    Triplet_Allocate(mass, *nDOF, 10.0 * non_zero_allocated);

    *u = calloc(*nDOF, sizeof(double));
    *v = calloc(*nDOF, sizeof(double));
    *a = calloc(*nDOF, sizeof(double));
    *F = calloc(*nDOF, sizeof(double));
    *dirichletValues = calloc(*nDOF, sizeof(double));

    for (int i = 0; i < M->NPoints; i++) {
        if (M->PointMark[i] == DIRICHLET_BC) {
            (*dirichletValues)[3 * i] = 0.0;
            (*dirichletValues)[3 * i + 1] = 0.0;
            (*dirichletValues)[3 * i + 2] = 0.0;
        }
    }

    for (int i = 0; i < M->NPoints; i++) {
        int dof_index = 3 * i;
        double x = M->x[i];
        (*v)[dof_index] = 0.0;
        (*v)[dof_index + 1] = 0.0;
        (*v)[dof_index + 2] = 0.0;
    }

    // for (int i = 0; i < M->NPoints; i++) {
    //     int dof_index = 3 * i;
    //     // double x = M->x[i];
    //     (*F)[dof_index] = 0.0;
    //     (*F)[dof_index + 1] = 0.0;//- 2*x;
    //     (*F)[dof_index + 2] = 0.0;//- 2*x;
    // }

    Mesh_Discretize_Mass_Matrix(M, mass);    
    CSR_InitFromTriplet(M_csr, mass);

    memset(K_csr, 0, sizeof(csr));
    memset(D_csr, 0, sizeof(csr));
    
    double *F_ext = *F;
    printf("Test\n");
    NewmarkIntegrate_Nonlinear(M, M_csr, D_csr, F_ext, *u, *v, *a, dt, nSteps, maxNewton, tolerance, useDamping, saveVibration, saveData, config, psi_force_3D, *dirichletValues, F_magnitude, azimuth_deg, elevation_deg );
}

void CleanupNonlinearElasticityProblem_Dynamics(mesh *M, triplet *mass, triplet *stiff, triplet *dump, csr *M_csr, csr *K_csr, csr *D_csr, double *u, double *v, double *a, double *F, double *dirichletValues) {
    Triplet_Free(mass);

    CSR_Free(M_csr);
    CSR_Free(K_csr);
    CSR_Free(D_csr);

    free(u);
    free(v);
    free(a);
    free(F);
    free(dirichletValues);

    Mesh_Free(M);
}

