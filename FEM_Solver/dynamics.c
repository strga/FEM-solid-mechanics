#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "csr.h"
#include "mesh.h"
#include "gelement3D.h"
#include "saving.h"
#include "dynamics.h"
#include "Mesh_GetBndrSide.h"
#include "triplet.h"
#include "Nonlinear.h"

// #define ENERGY_NONLINEAR

void initializeSimulationConfig(SimulationConfig *config) {
    // snprintf(config->displacementFile, sizeof(config->displacementFile), "Vibration_data/Linear_point_displacement_3D_cl021_11021e_damped.dat"); //Displacemen / vibration of chosen point
    // snprintf(config->energyFile, sizeof(config->energyFile), "Vibration_data/Linear_energy_3D_cl021_11021e_damped.txt"); // Conservation of energy file
    // snprintf(config->energyDampingFile, sizeof(config->energyDampingFile), "Vibration_data/Linear_energy_3D_cl021_11021e_damped.txt"); // Conservaiton of damped system
    // snprintf(config->timeSeriesFile, sizeof(config->timeSeriesFile), "Vibration_data/Wing_FFT_Dynamics.dat"); //FFT required data set

    // snprintf(config->displacementFile, sizeof(config->displacementFile), "Wing_matrix_etc/Wing_displacement_xyz_0(-90)0.dat"); //Displacemen / vibration of chosen point
    // snprintf(config->energyFile, sizeof(config->energyFile), "Wing_matrix_etc/Wing_conserve_energy_xyz_0(-90)0.txt"); // Conservation of energy file
    // snprintf(config->energyDampingFile, sizeof(config->energyDampingFile), "Wing_matrix_etc/Wing_conserve_damp_energy_xyz_0(-90)0.txt"); // Conservaiton of damped system
    // snprintf(config->timeSeriesFile, sizeof(config->timeSeriesFile), "Wing_matrix_etc/Wing_FFT_Dynamics_xyz_0(-90)0.dat"); //FFT required data set

    // FOR BEAMS 
    // snprintf(config->displacementFile, sizeof(config->displacementFile), "Wing_matrix_etc/Beam_displacement.dat"); //Displacemen / vibration of chosen point
    // snprintf(config->energyFile, sizeof(config->energyFile), "Wing_matrix_etc/Beam_energy_conserved.txt"); // Conservation of energy file
    // snprintf(config->energyDampingFile, sizeof(config->energyDampingFile), "Wing_matrix_etc/Beam_energy_conserved_damp.txt"); // Conservaiton of damped system
    // snprintf(config->timeSeriesFile, sizeof(config->timeSeriesFile), "Wing_matrix_etc/Beam_FFT.dat"); //FFT required data set

    // FOR WING
    snprintf(config->displacementFile, sizeof(config->displacementFile), "Beam_diplom4_displacement_damping_100K_small_new.dat"); //Displacemen / vibration of chosen point
    snprintf(config->energyFile, sizeof(config->energyFile), "Beam_diplom4_energy_conserved_100K_small_new.txt"); // Conservation of energy file
    snprintf(config->energyDampingFile, sizeof(config->energyDampingFile), "Beam_diplom4_energy_conserved_damping_100K_small_new.txt"); // Conservaiton of damped system
    snprintf(config->timeSeriesFile, sizeof(config->timeSeriesFile), "Beam_diplom4_FFT_damping_100K_small_new.dat"); //FFT required data set
}

int findNeumannBoundaryNode(const mesh *M) {
    int foundNode = -1;
    if (M->dim == 2) {
        for (int i = 0; i < M->NBEdges; i++) {
            Mesh_GetBndrSide S = MeshGetBoundarySide(M, i);

            // Check if S.idxA is on Neumann boundary
            if (M->PointMark[S.idxA] == NEUMANN_LOAD) {
                if (foundNode == -1 || S.idxA < foundNode) {
                    foundNode = S.idxA;
                }
            }

            // Check if S.idxB is on Neumann boundary
            if (M->PointMark[S.idxB] == NEUMANN_LOAD) {
                if (foundNode == -1 || S.idxB < foundNode) {
                    foundNode = S.idxB;
                }
            }
        }
    } else if (M->dim == 3) {
        for (int elem = 0; elem < M->NTriangles + M->NQuads; elem++) {
            if (M->ElementMark[elem] == NEUMANN_LOAD) {

                gelement3D face;
                GetElement3D(M, &face, elem);
                int nVertices = GetNumBasisFunctions(face.type);

                for (int v = 0; v < nVertices; v++) {
                    int node = face.idxNode[v];
                    if (foundNode == -1 || node < foundNode) {
                        foundNode = node;
                    }
                }
            }
        }
    } else {
        printf("Error: Unsupported mesh dimension (M->dim = %d)\n", M->dim);
        exit(EXIT_FAILURE);
    }
    
    return foundNode;
}

double calculateWorkByReference(double *u, double *U, csr *K, int nDOF) {
    double work_ref = 0.0, temp;
    for (int i = 0; i < nDOF; i++) {
        temp = 0.0;
        for (int j = K->PI[i]; j < K->PI[i + 1]; j++) {
            int colIndex = K->J[j];
            temp += K->VAL[j] * U[colIndex]; //U == reference displacement
        }
        work_ref += u[i] * temp; // result by current displacement 
    }
    return work_ref;
}

double calculateInitialEnergy(double *u_s, csr *K, int nDOF) {
    double initial_energy = 0.0, temp;
    for (int i = 0; i < nDOF; i++) {
        temp = 0.0;
        for (int j = K->PI[i]; j < K->PI[i + 1]; j++) {
            int colIndex = K->J[j];
            temp += K->VAL[j] * u_s[colIndex];
        }
        initial_energy += 0.5 * u_s[i] * temp;
    }
    return initial_energy;
}

double calculateKineticEnergy(double *v, csr *M, int nDOF) {
    double kinetic_energy = 0.0, temp;
    for (int i = 0; i < nDOF; i++) {
        temp = 0.0;
        for (int j = M->PI[i]; j < M->PI[i + 1]; j++) {
            int colIndex = M->J[j];
            temp += M->VAL[j] * v[colIndex];
        }
        kinetic_energy += 0.5 * v[i] * temp;
    }
    return kinetic_energy;
}

double calculatePotentialEnergy(double *u, csr *K, int nDOF) {
    double potential_energy = 0.0, temp;
    for (int i = 0; i < nDOF; i++) {
        temp = 0.0;
        for (int j = K->PI[i]; j < K->PI[i + 1]; j++) {
            int colIndex = K->J[j];
            temp += K->VAL[j] * u[colIndex];
        }
        potential_energy += 0.5 * u[i] * temp;
    }
    return potential_energy;
}

void NewmarkIntegrate(mesh *this, csr *M_csr, csr *K_csr, csr *D_csr , double *F, double *u, double *v, double *a, double dt, int nSteps, int useDamping, int saveVibration, int saveData, SimulationConfig *config) {
    
    int New_damp_energy_computaiton = TRUE;

    initializeSimulationConfig(config);
    const int nDOF = this->NPoints * this->dim;

    // Constants for Newmark integration
    const double beta = 0.25, gamma = 0.5;
    int neumannNode = findNeumannBoundaryNode(this);
    
    int* isDirichletNode = NULL;
    if (this->dim == 3) {
        isDirichletNode = calloc(this->NPoints, sizeof(int)); // 1 if node is constrained, 0 otherwise
        for (int elem = 0; elem < this->NTriangles + this->NQuads; elem++) {
            if (this->ElementMark[elem] == DIRICHLET_BC) {
                gelement3D K;
                GetElement3D(this, &K, elem);
                for (int i = 0; i < K.Nvertex; i++) {
                    int node = K.idxNode[i];
                    isDirichletNode[node] = 1;
                }
            }
        }
    }

    // Initialize CSR structures for damping and effective system matrix 
    csr A_csr;
    memset(&A_csr, 0, sizeof(csr));

    csr K_scaled;
    CSR_ScalarMultiply(K_csr, beta * dt * dt, &K_scaled); // beta * dt^2 * K
    CSR_Add(M_csr, &K_scaled, &A_csr);               // A = M + beta * dt^2 * K
    CSR_Free(&K_scaled);

    if (useDamping) {
        csr D_scaled;
        CSR_ScalarMultiply(D_csr, gamma * dt, &D_scaled); // gamma * dt * D

        csr temp_A;
        CSR_Add(&A_csr, &D_scaled, &temp_A);

        // Swap contents of temp_A with A_csr and free temp_A
        CSR_Free(&A_csr);
        A_csr = temp_A;

        CSR_Free(&D_scaled);
    }

    // Step 3: Allocate temporary vectors
    double *tempU = calloc(nDOF, sizeof(double));
    double *tempV = calloc(nDOF, sizeof(double));
    double *tempA = calloc(nDOF, sizeof(double));
    double *effectiveF = calloc(nDOF, sizeof(double));
    double *dampingForce = useDamping ? calloc(nDOF, sizeof(double)) : NULL;

    // Open files for saving data
    FILE *time_series_file = fopen(config->timeSeriesFile, "w");
    FILE *displacement_file = fopen(config->displacementFile, "w");
    FILE *energy_file = fopen(config->energyFile, "w");
    FILE *energy_dumping_file = fopen(config->energyDampingFile, "w");
    if (!time_series_file || !displacement_file || !energy_file || !energy_dumping_file) {
        fprintf(stderr, "Failed to open output files\n");
        if (time_series_file) fclose(time_series_file);
        if (displacement_file) fclose(displacement_file);
        if (energy_file) fclose(energy_file);
        if (energy_dumping_file) fclose(energy_dumping_file);
        return;
    }

    // Step 4: Initialize acceleration (a = M^-1 * (F - K * u))n
    CSR_Multiply(K_csr, u, effectiveF); // effectiveF = K * u
    if (useDamping) {
        CSR_Multiply(D_csr, v, dampingForce); // Initialize damping force
        for (int i = 0; i < nDOF; i++) {
            effectiveF[i] -= dampingForce[i];
        }
    }    
    for (int i = 0; i < nDOF; i++) {
        effectiveF[i] = F[i] - effectiveF[i];
    }
    CSR_Solve(M_csr, effectiveF, a);
    
    // Precompute initial energy correction
    double initial_energy = calculateInitialEnergy(u, K_csr, nDOF);
    
    double work_cumulative = 0.0; // ADDED

    // Newmark time integration loop
    for (int step = 0; step < nSteps; step++) {
        
        if (step % 1000 == 0) {
            printf("Done %d iterations\n", step);
        }    

        // Predictor step
        for (int i = 0; i < nDOF; i++) {
            tempU[i] = u[i] + dt * v[i] + 0.5 * dt * dt * (1 - 2 * beta) * a[i];
            tempV[i] = v[i] + dt * (1 - gamma) * a[i];
        }

        // Compute effective force
        CSR_Multiply(K_csr, tempU, effectiveF); // effectiveF = K * tempU
        for (int i = 0; i < nDOF; i++) {
            effectiveF[i] = F[i] - effectiveF[i];
        }

        if (useDamping) {
            CSR_Multiply(D_csr, tempV, dampingForce); // dampingForce = D * tempV
            for (int i = 0; i < nDOF; i++) {
                effectiveF[i] -= dampingForce[i];
            }
        }

        // Solve for acceleration
        CSR_Solve(&A_csr, effectiveF, tempA);

        // Corrector step
        for (int i = 0; i < nDOF; i++) {
            u[i] = tempU[i] + beta * dt * dt * tempA[i];
            v[i] = tempV[i] + gamma * dt * tempA[i];
            a[i] = tempA[i];
        }

        if (this->dim == 3 && isDirichletNode) {
            for (int node = 0; node < this->NPoints; ++node) {
                if (isDirichletNode[node]) {
                    for (int d = 0; d < this->dim; ++d) {
                        int dof = this->dim * node + d;
                        u[dof] = 0.0;
                        v[dof] = 0.0;
                        a[dof] = 0.0;
                    }
                }
            }
        }
        //----------------------------------// ADDED
        if (New_damp_energy_computaiton == TRUE) {
            double dotFv = 0.0;
            for (int i = 0; i < nDOF; i++) {
                dotFv += F[i] * v[i];
            }
            work_cumulative += dotFv * dt;
        }
        //----------------------------------//

        if (saveVibration && step < 80000) {
            char filename[256];
            snprintf(filename, sizeof(filename), "Newmark/100K_damping/Bend_damping100K_small_%06d.vtk", step); // Animaition of vibraiotn simulaton
            // snprintf(filename, sizeof(filename), "Wing_testing_2/Wing_-45_%06d.vtk", step); // Animaition of vibraiotn simulaton wing
            // snprintf(filename, sizeof(filename), "Wing_animation/Beam_silicon_%06d.vtk", step); // Animaition of vibraiotn simulaton
            
            if (this->dim == 2) {
                exportDisplacementToVTK(this, u, filename);
            } else if (this->dim == 3) {
                exportDisplacementToVTK_3D(this, u, filename);
            } else {
                printf("Error: Unsupported mesh dimension (this->dim = %d)\n", this->dim);
                exit(EXIT_FAILURE);
            }
        }

        // Save results
        if (saveData) {
            double time = (double)step * dt;
            // fprintf(displacement_file, "%.5f, %.5f, %.5f\n", time, u[2 * neumannNode], u[ 2 * neumannNode + 1]); // Displacement of point

            fprintf(displacement_file, "%.5f", time);  // Displacement of point
            for (int d = 0; d < this->dim; ++d)
                fprintf(displacement_file, ", %.5f", u[neumannNode * this->dim + d]);
            fprintf(displacement_file, "\n");

            double kinetic_energy = calculateKineticEnergy(v, M_csr, nDOF);
            double potential_energy = calculatePotentialEnergy(u, K_csr, nDOF);
            // double work_ref = calculateWorkByReference(u, F, K_csr, nDOF);

            if (!useDamping) {
                fprintf(energy_file, "%.5f, %.5f, %.5f, %.5f\n", time, kinetic_energy, potential_energy, kinetic_energy + potential_energy); // Conservation of energy
            }
            // if (useDamping) {
            //     fprintf(energy_dumping_file, "%.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n", time, kinetic_energy, potential_energy, initial_energy, work_ref, kinetic_energy + potential_energy - work_ref + initial_energy); // Conservation of energy
            // }
            if (useDamping) {
                fprintf(energy_dumping_file, "%.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n", time, kinetic_energy, potential_energy, initial_energy, work_cumulative, kinetic_energy + potential_energy - work_cumulative + initial_energy); // Conservation of energy ALTERNATIVE
            }
            if (step % 2 == 0) {
                if (this->dim == 2) {
                    saveTimeSeriesData(time_series_file, time, u, v, a, nDOF); // Data for FFT
                } else if (this->dim == 3) {
                    saveTimeSeriesData3D(time_series_file, time, u, v, a, nDOF); // Data for FFT
                } else {
                    printf("Error: Unsupported mesh dimension (this->dim = %d)\n", this->dim);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    // Close files
    fclose(displacement_file);
    fclose(energy_file);
    fclose(energy_dumping_file);
    fclose(time_series_file);

    CSR_Free(&A_csr);
    free(tempU); free(tempV); free(tempA); free(effectiveF);
    if (isDirichletNode) free(isDirichletNode);
    if (dampingForce) free(dampingForce);
}

void NewmarkIntegrate_Nonlinear(mesh *M, csr *M_csr, csr *D_csr, double *F_ext, double *u, double *v, double *a, double dt, int nSteps, int maxNewton, double tol, int useDamping, int saveVibration, int saveData, SimulationConfig *config, double (*psi_force_3D)(double,double,double,double,double,double,double,double,double,double*,double*,double*), double *dirichletValues, double F_mag, double azimuth_deg, double elevation_deg)
{
    const int nDOF = M->dim * M->NPoints;
    const double beta = 0.25, gamma = 0.5;
    const double a0 = 1.0 / (beta * dt * dt);
    const double a1 = gamma / (beta * dt);

    // Allocate persistent work arrays
    double *u_pred = calloc(nDOF, sizeof(double));
    double *v_pred = calloc(nDOF, sizeof(double));
    double *a_trial = calloc(nDOF, sizeof(double));
    double *v_trial = calloc(nDOF, sizeof(double));
    double *delta_u = calloc(nDOF, sizeof(double));
    double *residual = calloc(nDOF, sizeof(double));
    double *P_int = calloc(nDOF, sizeof(double));
    double *Ma = calloc(nDOF, sizeof(double));
    double *Dv = useDamping ? calloc(nDOF, sizeof(double)) : NULL;

    triplet K_triplet;
    int nnz_guess = (int)(4.5 * M->NTetrahedra * 144); // adapt based on mesh
    Triplet_Init(&K_triplet, nDOF, nDOF, nnz_guess);

    csr K_csr = {0}, A_csr = {0};

    for (int step = 0; step < nSteps; step++) {
        if (step % 10 == 0)
            printf("Time step %d\n", step);

        // Predictor step
        for (int i = 0; i < nDOF; i++) {
            u_pred[i] = u[i] + dt * v[i] + 0.5 * dt * dt * (1 - 2 * beta) * a[i];
            v_pred[i] = v[i] + dt * (1 - gamma) * a[i];
        }

        // Initialize u to predictor for Newton iterations
        memcpy(u, u_pred, nDOF * sizeof(double));

        // Newton-Raphson loop
        for (int iter = 0; iter < maxNewton; iter++) {

            // Assemble internal force and tangent stiffness
            Triplet_Reset(&K_triplet);
            Mesh_Discretize_NonlinearDynamics3D(M, &K_triplet, residual, F_ext, P_int, u, dirichletValues, psi_force_3D, F_mag, azimuth_deg, elevation_deg);

            if (K_csr.nz > 0)
                CSR_Free(&K_csr);
            CSR_InitFromTriplet(&K_csr, &K_triplet);

            // Build effective tangent matrix: A = a0*M + a1*D + K
            csr M_scaled, D_scaled, A_temp;
            CSR_ScalarMultiply(M_csr, a0, &M_scaled);
            if (useDamping)
                CSR_ScalarMultiply(D_csr, a1, &D_scaled);

            CSR_Add(&M_scaled, &K_csr, &A_temp);
            if (useDamping) {
                csr A_temp2;
                CSR_Add(&A_temp, &D_scaled, &A_temp2);
                CSR_Free(&A_temp);
                A_temp = A_temp2;
                CSR_Free(&D_scaled);
            }
            CSR_Free(&M_scaled);

            if (A_csr.nz > 0)
                CSR_Free(&A_csr);
            A_csr = A_temp;

            // Compute trial acceleration and velocity based on current u
            for (int i = 0; i < nDOF; i++) {
                a_trial[i] = a0 * (u[i] - u_pred[i]);
                v_trial[i] = v_pred[i] + dt * gamma * a_trial[i];
            }

            // Compute effective residual
            CSR_Multiply(M_csr, a_trial, Ma);
            if (useDamping)
                CSR_Multiply(D_csr, v_trial, Dv);

            for (int i = 0; i < nDOF; i++) {
                residual[i] = Ma[i] + (useDamping ? Dv[i] : 0.0) + P_int[i] - F_ext[i];
            }

            // Solve linear system for correction delta_u
            CSR_Solve(&A_csr, residual, delta_u);

            // Update u
            for (int i = 0; i < nDOF; i++)
                u[i] -= delta_u[i];

            // Check convergence
            double norm_du = 0.0, norm_R = 0.0;
            for (int i = 0; i < nDOF; i++) {
                norm_du += delta_u[i] * delta_u[i];
                norm_R += residual[i] * residual[i];
            }
            norm_du = sqrt(norm_du);
            norm_R = sqrt(norm_R);

            printf("   Newton iter %d: ||du||=%.3e, ||R||=%.3e\n", iter, norm_du, norm_R);
            if (norm_du < tol && norm_R < tol)
                break;
        }

        // Corrector step: update acceleration and velocity after Newton convergence
        for (int i = 0; i < nDOF; i++) {
            a[i] = a0 * (u[i] - u_pred[i]);
            v[i] = v_pred[i] + dt * gamma * a[i];
        }

        // Optional export
        if (saveVibration) {
            char filename[256];
            printf("VTK saved\n");
            snprintf(filename, sizeof(filename), "nonlinear_step_%05d.vtk", step);
            exportDisplacementToVTK_3D(M, u, filename);
        }
    }

    // Cleanup
    Triplet_Free(&K_triplet);
    CSR_Free(&K_csr);
    CSR_Free(&A_csr);
    free(u_pred); free(v_pred); free(a_trial); free(v_trial);
    free(delta_u); free(residual); free(P_int); free(Ma);
    if (Dv) free(Dv);
}
