/* ----------------------------------------------------------------------------------------- */
#include <stdio.h>
#include "constants.h"
/* ----------------------------------------------------------------------------------------- */
#ifndef _GMSHLINE_H
#define _GMSHLINE_H
/* ----------------------------------------------------------------------------------------- */
#define MESHFORMAT_MSH  2
#define MESHFORMAT_MSH4 4
#define MAX_PHYSICAL   20
#define MAX_SUBDOMAINS 500
typedef struct {
    int idx;
    int etyp;
    int etags;
    int tagPhysical;
    int tagGeometrical;
    int tagSubDomain;
    int ilist[50]; 
} gmshline;	

typedef struct {
    int NPoints;
    int NGeomElements;
    int NBndrEdges;
    int NElements;
    int NTriangles;
    int NQuads;
    int NTetrahedra;  // Number of tetrahedra
    int NHexahedra;   // Number of hexahedra
    int NWedges;      // Number of wedges (prisms)
    int NPyramids;  // Number of pyramids
    int NInterfaceEdges;
    int NEdges;
    /* Physical domains */
    int nPhysical;
    int labelPhysical[MAX_PHYSICAL];
    int ntriPhysical[MAX_PHYSICAL];
    int nquadPhysical[MAX_PHYSICAL];
    int nedgePhysical[MAX_PHYSICAL];
    int ntetPhysical[MAX_PHYSICAL];  // Tetrahedra
    int nhexPhysical[MAX_PHYSICAL];  // Hexahedra
    int nwedgePhysical[MAX_PHYSICAL]; // Wedges
    int npyrPhysical[MAX_PHYSICAL];  // Pyramids
    /* Subdomains */
    int nSubdomain;
    int labelSubdomain[MAX_SUBDOMAINS];
    int ntriSubdomain[MAX_SUBDOMAINS];
    int nquadSubdomain[MAX_SUBDOMAINS];
    int nedgeSubdomain[MAX_SUBDOMAINS];
    int ntetSubdomain[MAX_SUBDOMAINS];  // Tetrahedra
    int nhexSubdomain[MAX_SUBDOMAINS];  // Hexahedra
    int nwedgeSubdomain[MAX_SUBDOMAINS]; // Wedges
    int npyrSubdomain[MAX_SUBDOMAINS]; // Pyramid
} gmshfile;
/* ----------------------------------------------------------------------------------------- */
int        gmshfile_detectformat(const char *fname);
/* ----------------------------------------------------------------------------------------- */
void       gmshfile_read(gmshfile *p, const char *fname);
int        gmshline_read(gmshline *p, const char *buf);
void       gmshfile_add2subdomain(gmshfile *p, int id, int etyp);
void       gmshfile_add2physical(gmshfile *p, int mark, int etyp);
int        gmsh_getNPoints(const char *fname);
int        gmsh_getNodesCoordinates(FILE *fid, double *x, double *y, double *z, int dim);
/* ------------------------------------------------------------------------------------------ */
#endif
















