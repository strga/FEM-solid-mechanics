#include <stdio.h>
#include <stdlib.h>
#include "saving.h"
#include "gelement.h"
#include "gelement3D.h"

void WriteTripletMatrixToFile(const triplet *p, const char *filename) {
    FILE *fid = fopen(filename, "w+");
    if (fid == NULL) {
        fprintf(stderr, "Error opening file to write matrix\n");
        return;
    }

    for (int i = 0; i < p->nz; i++) {
        fprintf(fid, "%10d %10d %24.15g\n", p->tripletI[i] + 1, p->J[i] + 1, p->VAL[i]);
    }

    fclose(fid);
}

void exportNumericalSolutionToVTK(const mesh *M, const double *numericalSolution, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Numerical solution\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(file, "POINTS %d float\n", M->NPoints);
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f %f %f\n", M->x[i], M->y[i], numericalSolution[i]); 
    }

    fprintf(file, "CELLS %d %d\n", M->NElements, 4 * M->NElements);
    for (int i = 0; i < M->NElements; i++) {
        gelement K = Gelement(M, i);
        fprintf(file, "3 %d %d %d\n", K.idx[0], K.idx[1], K.idx[2]);
    }

    fprintf(file, "CELL_TYPES %d\n", M->NElements);
    for (int i = 0; i < M->NElements; i++) {
        fprintf(file, "5\n");
    }

    fprintf(file, "POINT_DATA %d\n", M->NPoints);
    fprintf(file, "SCALARS solution float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f\n", numericalSolution[i]);
    }

    fclose(file);
}

void exportNumericalSolutionToVTK_3D(const mesh *M, const double *numericalSolution, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Numerical solution\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    // Write 3D points (x, y, z) from the mesh
    fprintf(file, "POINTS %d float\n", M->NPoints);
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f %f %f\n", M->x[i], M->y[i], M->z[i]);
    }

    // Index ranges for each element type
    int startQuad = M->NTriangles;
    int startTetra = startQuad + M->NQuads;
    int startHex = startTetra + M->NTetrahedra;
    int startPrism = startHex + M->NHexahedra;
    int startPyramid = startPrism + M->NWedges;

    // Count the 3D elements
    int nVolumeElements = M->NTetrahedra + M->NHexahedra + M->NPyramids + M->NWedges;
    int nIndex = M->NTetrahedra * 5 + M->NHexahedra * 9 + M->NPyramids * 6 + M->NWedges * 7; 
    // Write cells with the correct node counts for 3D elements
    int cell_count = 0;
    fprintf(file, "CELLS %d %d\n", nVolumeElements, nIndex);
    for (int i = M->NTriangles + M->NQuads; i < M->NElements; i++) {
        gelement3D K; 
        GetElement3D(M, &K, i);

        if (i >= startTetra && i < startHex && K.Nvertex == 4) {
            // Tetrahedron (4 vertices)
            fprintf(file, "4 %d %d %d %d\n",
                    K.idxNode[0], K.idxNode[1], K.idxNode[2], K.idxNode[3]);
            cell_count++;
        } 
        else if (i >= startHex && i < startPrism && K.Nvertex == 8) {
            // Hexahedron (8 vertices)
            fprintf(file, "8 %d %d %d %d %d %d %d %d\n",
                    K.idxNode[0], K.idxNode[1], K.idxNode[2], K.idxNode[3],
                    K.idxNode[4], K.idxNode[5], K.idxNode[6], K.idxNode[7]);
            cell_count++;
        } 
        else if (i >= startPrism && i < startPyramid && K.Nvertex == 6) {
            // Wedge/Prism (6 vertices)
            fprintf(file, "6 %d %d %d %d %d %d\n",
                    K.idxNode[0], K.idxNode[1], K.idxNode[2],
                    K.idxNode[3], K.idxNode[4], K.idxNode[5]);
            cell_count++;
        } 
        else if (i >= startPyramid && K.Nvertex == 5) {
            // Pyramid (5 vertices)
            fprintf(file, "5 %d %d %d %d %d\n",
                    K.idxNode[0], K.idxNode[1], K.idxNode[2], K.idxNode[3],
                    K.idxNode[4]);
            cell_count++;
        }
    }

    // Write cell types
    fprintf(file, "CELL_TYPES %d\n", cell_count);
    for (int i = startTetra; i < M->NElements; i++) {
        if (i >= startTetra && i < startHex) {
            fprintf(file, "10\n"); // VTK_TETRAHEDRON
        } 
        else if (i >= startHex && i < startPrism) {
            fprintf(file, "12\n"); // VTK_HEXAHEDRON
        } 
        else if (i >= startPrism && i < startPyramid) {
            fprintf(file, "13\n"); // VTK_WEDGE (Prism)
        } 
        else if (i >= startPyramid) {
            fprintf(file, "14\n"); // VTK_PYRAMID
        }
    }

    // Write scalar solution data for each point
    fprintf(file, "POINT_DATA %d\n", M->NPoints);
    fprintf(file, "SCALARS solution float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f\n", numericalSolution[i]);
    }

    fclose(file);
}

void exportDisplacementToVTK(const mesh *M, const double *displacement, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\nDisplacement Visualization\nASCII\nDATASET UNSTRUCTURED_GRID\n");

    fprintf(file, "POINTS %d float\n", M->NPoints);
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f %f 0.0\n", M->x[i], M->y[i]); 
    }

    fprintf(file, "CELLS %d %d\n", M->NElements, 4 * M->NElements);
    for (int i = 0; i < M->NElements; i++) {
        gelement K = Gelement(M, i);
        fprintf(file, "3 %d %d %d\n", K.idx[0], K.idx[1], K.idx[2]);
    }

    fprintf(file, "CELL_TYPES %d\n", M->NElements);
    for (int i = 0; i < M->NElements; i++) {
        fprintf(file, "5\n");
    }

    fprintf(file, "POINT_DATA %d\n", M->NPoints);
    fprintf(file, "VECTORS displacement_vectors float\n");
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f %f 0.0\n", displacement[2*i], displacement[2*i+1]);
    }

    fclose(file);
}

void exportDisplacementToVTK_3D(const mesh *M, const double *displacement, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "3D Vector Displacement\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    // Write mesh points
    fprintf(file, "POINTS %d float\n", M->NPoints);
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f %f %f\n", M->x[i], M->y[i], M->z[i]);
    }

    // Write 3D elements (similar to your existing scalar exporter)
    int startQuad = M->NTriangles;
    int startTetra = startQuad + M->NQuads;
    int startHex = startTetra + M->NTetrahedra;
    int startPrism = startHex + M->NHexahedra;
    int startPyramid = startPrism + M->NWedges;

    int nVolumeElements = M->NTetrahedra + M->NHexahedra + M->NPyramids + M->NWedges;
    int nIndex = M->NTetrahedra * 5 + M->NHexahedra * 9 + M->NPyramids * 6 + M->NWedges * 7;

    fprintf(file, "CELLS %d %d\n", nVolumeElements, nIndex);
    for (int i = M->NTriangles + M->NQuads; i < M->NElements; i++) {
        gelement3D K; 
        GetElement3D(M, &K, i);

        if (i >= startTetra && i < startHex && K.Nvertex == 4) {
            fprintf(file, "4 %d %d %d %d\n", K.idxNode[0], K.idxNode[1], K.idxNode[2], K.idxNode[3]);
        } 
        else if (i >= startHex && i < startPrism && K.Nvertex == 8) {
            fprintf(file, "8 %d %d %d %d %d %d %d %d\n",
                K.idxNode[0], K.idxNode[1], K.idxNode[2], K.idxNode[3],
                K.idxNode[4], K.idxNode[5], K.idxNode[6], K.idxNode[7]);
        } 
        else if (i >= startPrism && i < startPyramid && K.Nvertex == 6) {
            fprintf(file, "6 %d %d %d %d %d %d\n",
                K.idxNode[0], K.idxNode[1], K.idxNode[2],
                K.idxNode[3], K.idxNode[4], K.idxNode[5]);
        } 
        else if (i >= startPyramid) {
            fprintf(file, "5 %d %d %d %d %d\n",
                K.idxNode[0], K.idxNode[1], K.idxNode[2], K.idxNode[3],
                K.idxNode[4]);
        }
    }

    // Write cell types
    fprintf(file, "CELL_TYPES %d\n", nVolumeElements);
    for (int i = M->NTriangles + M->NQuads; i < M->NElements; i++) {
        if (i >= startTetra && i < startHex) {
            fprintf(file, "10\n"); // VTK_TETRAHEDRON
        } 
        else if (i >= startHex && i < startPrism) {
            fprintf(file, "12\n"); // VTK_HEXAHEDRON
        } 
        else if (i >= startPrism && i < startPyramid) {
            fprintf(file, "13\n"); // VTK_WEDGE
        } 
        else if (i >= startPyramid) {
            fprintf(file, "14\n"); // VTK_PYRAMID
        }
    }

    // Write 3D displacement vectors
    fprintf(file, "POINT_DATA %d\n", M->NPoints);
    fprintf(file, "VECTORS displacement float\n");
    for (int i = 0; i < M->NPoints; i++) {
        fprintf(file, "%f %f %f\n", 
            displacement[3 * i], 
            displacement[3 * i + 1], 
            displacement[3 * i + 2]);
    }

    fclose(file);
}

void saveVectorToFile(double *vector, int size, int dimensions, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Failed to open file %s for writing.\n", filename);
        return;
    }

    if (dimensions == 1) {
        for (int i = 0; i < size; i++) {
            fprintf(file, "%.15f\n", vector[i]);
        }
    } else if (dimensions == 2) {
        for (int i = 0; i < size; i += 2) {
            fprintf(file, "%.15e %.15e\n", vector[i], vector[i + 1]);
        }
    } else if (dimensions == 3) {
        for (int i = 0; i < size; i += 3) {
            fprintf(file, "%.15e %.15e %.15e\n", vector[i], vector[i + 1], vector[i + 2]);
        }
    }

    fclose(file);
}

void saveSparseMatrixToFile(const csr *A_csr, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Failed to open file %s for writing.\n", filename);
        return;
    }

    for (int i = 0; i < A_csr->n; i++) {
        for (int j = A_csr->PI[i]; j < A_csr->PI[i + 1]; j++) {
            fprintf(file, "%d %d %.15e\n", i + 1, A_csr->J[j] + 1, A_csr->VAL[j]); // +1 for MATLAB compatibility
        }
    }

    fclose(file);
}

void saveTimeStepData(const char *filename, double *displacement, double *velocity, int NPoints, double time) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "Time: %.5f\n", time);
    for (int i = 0; i < NPoints; i++) {
        fprintf(file, "%.15e %.15e %.15e %.15e\n", displacement[2 * i], displacement[2 * i + 1], velocity[2 * i], velocity[2 * i + 1]);
    }

    fclose(file);
}

void saveTimeSeriesData(FILE *time_series_file, double time, double *u, double *v, double *a, int nDOF) {
    if (!time_series_file) {
        fprintf(stderr, "Error: Invalid file pointer for saving time-series data.\n");
        return;
    }

    fprintf(time_series_file, "%.15f", time);
    for (int i = 0; i < nDOF / 2; i++) {
        int x_index = 2 * i;
        int y_index = 2 * i + 1;

        // fprintf(time_series_file, ", %.15f, %.15f, %.15f, %.15f, %.15f, %.15f", u[x_index], u[y_index], v[x_index], v[y_index], a[x_index], a[y_index]);
        fprintf(time_series_file, ", %.15f, %.15f", u[x_index], u[y_index]);
    }

    fprintf(time_series_file, "\n");
}

void saveTimeStepData3D(const char *filename, double *displacement, double *velocity, int NPoints, double time) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "Time: %.5f\n", time);
    for (int i = 0; i < NPoints; i++) {
        fprintf(file, "%.15e %.15e %.15e %.15e %.15e %.15e\n", displacement[3 * i], displacement[3 * i + 1], displacement[3 * i + 2], velocity[3 * i], velocity[3 * i + 1], velocity[3 * i + 2]);
    }

    fclose(file);
}

void saveTimeSeriesData3D(FILE *time_series_file, double time, double *u, double *v, double *a, int nDOF) {
    if (!time_series_file) {
        fprintf(stderr, "Error: Invalid file pointer for saving time-series data.\n");
        return;
    }

    fprintf(time_series_file, "%.15f", time);
    for (int i = 0; i < nDOF / 3; i++) {
        int x_index = 3 * i;
        int y_index = 3 * i + 1;
        int z_index = 3 * i + 2;

        // fprintf(time_series_file, ", %.15f, %.15f, %.15f, %.15f, %.15f, %.15f", u[x_index], u[y_index], v[x_index], v[y_index], a[x_index], a[y_index]);
        fprintf(time_series_file, ", %.15f, %.15f, %.15f", u[x_index], u[y_index], u[z_index]);
    }

    fprintf(time_series_file, "\n");
}

void saveIterationValue(const char *filename, int iteration, double value) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    fprintf(file, "%d %.15e\n", iteration, value);
    fclose(file);
}

