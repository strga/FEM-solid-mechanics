# FEM Solid Mechanics Solver

3D finite element solver for solid mechanics problems, implemented in C. Supports linear and nonlinear static and dynamic analysis for isotropic and orthotropic materials, with MATLAB post-processing scripts.

## Capabilities

| Mode | Description |
|------|-------------|
| `POISSON2D` / `POISSON3D` | Verification: Laplace/Poisson equation with analytical solution |
| `STATICS2D` / `STATICS3D` | Linear elasticity, static loading |
| `STATICS3D_ORTHO` | Linear elasticity, orthotropic materials |
| `DYNAMICS2D` / `DYNAMICS3D` | Newmark-β time integration, isotropic |
| `DYNAMICS3D_ORTHO` | Newmark-β time integration, orthotropic |
| `NONLINEAR` | Nonlinear statics — Saint Venant-Kirchhoff and Neo-Hookean hyperelasticity |
| `NONLINEAR_DYNAMICS` | Nonlinear dynamics with Newton iteration - UNFINISHED|

## Methods

- Galerkin FEM with isoparametric elements (triangles, quads, tetrahedra, hexahedra, wedges, pyramids)
- Sparse stiffness and mass matrix assembly in CSR format
- Direct solver: UMFPACK
- Iterative solver: PETSc KSP
- Newmark-beta scheme for time integration
- Newton-Raphson iteration for nonlinear problems
- VTK output for visualization in ParaView

## Repository Structure
