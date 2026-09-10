/* ----------------------------------------------------------------------------------------- */
#ifndef _CONSTANTS_H
#define _CONSTANTS_H
/* ----------------------------------------------------------------------------------------- */
#define VTK_LINE                3
#define VTK_TRIANGLE            5
#define VTK_QUAD                9
#define VTK_TETRA               10
#define VTK_HEXAHEDRON          12
#define VTK_WEDGE               13
#define VTK_PYRAMID             14
/* ----------------------------------------------------------------------------------------- */
#define GMSH_POINT            15
#define GMSH_SEGMENT          1
#define GMSH_TRIANGLE         2
#define GMSH_QUADRANGLE       3
#define GMSH_TETRAHEDRON      4
#define GMSH_HEXAHEDRON       5
#define GMSH_WEDGE            6
#define GMSH_PYRAMID          7
/* ----------------------------------------------------------------------------------------- */
#define SQU(a)               ((a)*(a))
#define max(a,b)             (((a) > (b)) ? (a) : (b) )
#define min(a,b)             (((a) < (b)) ? (a) : (b) ) 
#define pair1(a,b)           (((a) < (b)) ? (a) : (b) ) 
#define pair2(a,b)           (((a)<(b)) ? (b) : (a) ) 
#define EuclidianNorm(a,b,c) (sqrt((a) * (a) + (b) * (b) + (c) * (c))) 
#define EuclidianNorm2(a,b)  (sqrt((a) * (a) + (b) * (b))) 
#define IsNonZero(x)         ((x > 1e-20) ? 1 : ((x < -1e-20) ? 1 : 0) )
#define IsZero(x)            ((x > 1e-20) ? 0 : ((x < -1e-20) ? 0 : 1) )
/* ----------------------------------------------------------------------------------------- */
#define maxf(a,b,c,d)         (max(max(a, b), max(c, d)))
#define minf(a,b,c,d)         (min(min(a, b), min(c, d)))
/* ----------------------------------------------------------------------------------------- */
typedef enum { BDF1, BDF2, BDF3, BDF4 } timescheme;
/* ----------------------------------------------------------------------------------------- */
typedef enum {UNDEF_FE,          LINEAR_FE, QUADRATIC_FE, CUBIC_FE, P4_FE, P5_FE, P6_FE, P7_FE, P8_FE, P9_FE }      fetype;
typedef enum {UNDEF_FFE,         TAYLORHOOD_FFE, P3P2_FFE, P4P3_FFE, P1P1_FFE }                                     flowfetype;
typedef enum {UNDEF_FEBASETYPE,  FEBASETYPE_DEFAULT, FEBASETYPE_HIERARCHIC, FEBASETYPE_LAGRANGIAN }                 febasetype;
typedef enum {UNDEF_PROBLEM,  SCALAR_PROBLEM,     IFLOW_PROBLEM,       ELASTIC_PROBLEM}               problemtype;
/* ----------------------------------------------------------------------------------------- */
#define MAXIMUM_NUMBER_OF_MARKS 50
/* ----------------------------------------------------------------------------------------- */
#define MARK_DIRICHLET           1 
#define MARK_NEWTON             20
#define MARK_NEUMANN            60
/* ----------------------------------------------------------------------------------------- */
#define MARK_INLET               1 
#define MARK_MOVING_BOUNDARY    10
#define MARK_OUT                20
#define MARK_WALL               30
#define MARK_ZERO               40
#define MARK_SYMMETRY_XZplane   60
#define MARK_SLIP_BC_XZplane    81
#define MARK_SYMMETRY_YZplane   65
#define MARK_SLIP_BC_YZplane    82
#define MARK_PRESSURE_PEN       220

#define MARK_SYMMETRY_XYplane    63
#define MARK_SYMMETRY            70
#define MARK_SLIP_BC             80
#define MARK_SLIP_BC_XYplane     83
/* ----------------------------------------------------------------------------------------- */
#define MAXIMUM_TIMEHISTORY      3
#define MAX_BASEFUNCTIONS_2D     16
#define MAX_BASEFUNCTIONS_3D     16
/* ----------------------------------------------------------------------------------------- */
#define MAX_SYSTEM_DIM           5
#define MAX_LOCALSYSTEM          (MAX_BASEFUNCTIONS_2D)
#define MAX_VECTORLOCALSYSTEM    (MAX_SYSTEM_DIM * MAX_BASEFUNCTIONS_2D)
/* ----------------------------------------------------------------------------------------- */
#define BUFFER_SIZE              500 
#define FILENAME_LENGTH          128
/* ----------------------------------------------------------------------------------------- */
#define INDEX_NOT_SET -1
#define TRUE           1
#define FALSE          0
/* ----------------------------------------------------------------------------------------- */
#define MASTERELEMENT_REFINE             64
#define MAX_NCUT 200
/* ----------------------------------------------------------------------------------------- */
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
/* ----------------------------------------------------------------------------------------- */
#endif
