#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <umfpack.h>
#include "Nonlinear.h"
#include "csr.h"
#include "gelement3D.h"
#include "mesh.h"
#include "triplet.h"
#include "saving.h"

void runNonlinearWithLoadSteps(mesh *M, double *u_global, triplet *stiff, double *residual, double *f_ext, double *P_int, double *dirichletValues, int maxIter, double tolerance, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_total, double azimuth_deg, double elev_deg, int nSteps)
{
    for (int step = 1; step <= nSteps; ++step) {
        double F_step = F_total * ((double)step / nSteps);
        printf("\n========= Load Step %d / %d: Applied Load = %.4e =========\n", step, nSteps, F_step);

        solveNonLinear(M, u_global, stiff, residual, f_ext, P_int, dirichletValues, maxIter, tolerance, psi_force_3D, F_step, azimuth_deg, elev_deg);

        // Monitor deformation growth
        double u_max = 0.0;
        for (int i = 0; i < 3 * M->NPoints; ++i)
            u_max = fmax(u_max, fabs(u_global[i]));
        printf("→ Max displacement after step %d: %.6e\n", step, u_max);
    }
}

void solveNonLinear(mesh *M, double *u_global, triplet *stiff, double *residual, double *f_ext, double *P_int, double *dirichletValues, int maxIter, double tolerance, double (*psi_force_3D)(double, double, double, double, double, double, double, double, double, double*, double*, double*), double F_mag, double azimuth_deg, double elev_deg)
{
    const int nDOF = 3 * M->NPoints;
    double *du = calloc(nDOF, sizeof(double));
    double *rhs = calloc(nDOF, sizeof(double));
    csr K_csr;

    double norm_du0_sq = 1.0;

    FILE *residual_file = fopen("residual_R_log_case5.txt", "w");
    FILE *du_file = fopen("residual_du_log_case5.txt", "w");

    if (!residual_file || !du_file) {
        fprintf(stderr, "Failed to open residual output files.\n");
        if (residual_file) fclose(residual_file);
        if (du_file) fclose(du_file);
        return;
    }
    
    double prev_convergence_R = -1.0;
    const double residual_stagnation_tol = 1e-25;  // You can tune this value

    for (int iter = 0; iter < maxIter; iter++) {
        
        // if (iter % 100 == 0) { 
        //     printf("Iteration %d:\n", iter);
        // }
        printf("Iteration %d:\n", iter);
    
        Triplet_Reset(stiff);
        Mesh_Discretize_NonlinearStatics3D(M, stiff, residual, f_ext, P_int, u_global, dirichletValues, psi_force_3D, F_mag, azimuth_deg, elev_deg);
        CSR_InitFromTriplet(&K_csr, stiff);

        for (int i = 0; i < nDOF; i++) {
            rhs[i] = -residual[i];
        }
        CSR_Solve(&K_csr, rhs, du); // Solve system K * du = residual

        double norm_du_sq = 0.0;
        for (int i = 0; i < nDOF; i++) {
            norm_du_sq += du[i] * du[i];
        }
        if (iter == 0) {
            norm_du0_sq = norm_du_sq > 1e-25 ? norm_du_sq : 1.0;
        }
        double convergence_du = norm_du_sq / (1.0 + norm_du0_sq);


        for (int i = 0; i < nDOF; i++) {
            u_global[i] += du[i];
        }

        // Triplet_Reset(stiff);
        // Mesh_Discretize_NonlinearStatics3D(M, stiff, residual, f_ext, P_int, u_global, dirichletValues, psi_force_3D, F_mag, azimuth_deg, elev_deg);

        double RNumerator = 0.0;
        double RDenominator = 1.0;
        for (int i = 0; i < nDOF; i++) {
            RNumerator += residual[i] * residual[i];
            RDenominator += f_ext[i] * f_ext[i];
        }
        double convergence_R = sqrt( RNumerator / RDenominator );
        
        double delta_convergence_R = prev_convergence_R;


        fprintf(residual_file, "%d %.15e\n", iter, convergence_R);
        fprintf(du_file, "%d %.15e\n", iter, convergence_du);

        printf("  Convergence ||du||: %.15e\n", convergence_du);
        printf("  Convergence ||R||:  %.15e\n", convergence_R);
        printf("  Convergence ||dR||:  %.15e\n", delta_convergence_R);

        // if (delta_convergence_R < residual_stagnation_tol) {
        //     printf("Converged due to residual stagnation at iteration %d.\n", iter + 1);
        //     printf("  Residual change: %.6e\n", delta_convergence_R);
        //     break;
        // }

        // if (convergence_R < tolerance && convergence_du < tolerance) {
        //     printf("Converged after %d iterations.\n", iter + 1);
        //     printf("  Final convergence ||du||: %.6e\n", convergence_du);
        //     printf("  Final convergence ||R||:  %.6e\n", convergence_R);
        //     break;
        // }

        if (iter == maxIter - 1) {
            printf("  Newton-Raphson did not converge after %d iterations.\n", maxIter);
            printf("  Final ||R|| = %.6e, ||du|| = %.6e\n", convergence_R, convergence_du);
        }

        prev_convergence_R = convergence_R;

        // if (iter % 20 == 0) {
        //     printf("Convergence of R: %g\n", convergence_R);
        // }
        
        // if (convergence_R < tolerance) {
        //     printf("Convergence of R was satisfied: Final convergence is %g\n", convergence_R);
        //     break;
        // }

        // if (iter == maxIter - 1) {
        //     printf("Warning: Newton-Raphson did not converge after %d iterations.\n", maxIter);
        //     printf("Final residdum: %g\n", convergence_R);
        //     break;
        // }

        if (iter == 40)
        {
            break;
        }
        
        CSR_Free(&K_csr);
    }

    fclose(residual_file);
    fclose(du_file);

    free(du);
    free(rhs);
}

