/* ------------------------------------------------------------------------------------------ */
#include "constants.h"
#include "wlog.h"
#include "gmshtools.h"
/* ------------------------------------------------------------------------------------------ */
#include <string.h>
#include <strings.h>
#include <stdlib.h>
/* ------------------------------------------------------------------------------------------------------- */
int gmshfile_detectformat(const char *fname)
{
    /** @bug - to be corrected */
    return MESHFORMAT_MSH;
}
/* ------------------------------------------------------------------------------------------ */
int gmshline_read(gmshline *p, const char *buf)
{
    int idx, etp, tgs, tgs1;
    int val;
    int tags[10], pom[8];
    int nred;
    // printf("Raw buffer content: %s\n", buf);

    nred = sscanf(buf, "%d %d %d", &idx, &etp, &tgs);
    if (nred != 3)
    {
        writelog("Error[MSH]:\t GMSH format is wrong, wrong line or buffer!\n");
        exit(1);
    }
    p->etyp = etp;
    p->idx  = idx;
    tags[0] = tags[1] = tags[2] = tags[3] = tags[4]= tags[5] = pom[0] = pom[1] = pom[2] = pom[3] = 0;
    // printf("Parsed idx: %d, etp: %d, tgs: %d\n", idx, etp, tgs);

    p->tagSubDomain = INDEX_NOT_SET;

    val  = tgs + etp * 10;
    nred = 0;

    // printf("tgs1: %d \t nred1: %d\n",tgs, nred);

    switch (val)
    {
    case 11:
        nred=sscanf(buf, "%d %d %d %d %d %d"   , &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1]);
        break;
    case 12:
        nred=sscanf(buf, "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1]);
        break;
    case 13:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &pom[0], &pom[1]);
        break;
    case 14:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &pom[0], &pom[1]);
        break;
    case 15:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &pom[0], &pom[1]);
        break;
    case 16:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &tags[5], &pom[0], &pom[1]);
        break;
    case 17:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &tags[5], &tags[6], &pom[0], &pom[1]);
        break;
    case 21:
        nred=sscanf(buf, "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1], &pom[2]);
        break;
    case 22:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2]);
        break;
    case 23:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &pom[0], &pom[1], &pom[2]);
        break;
    case 24:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &pom[0], &pom[1], &pom[2]);
        break;
    case 25:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &pom[0], &pom[1], &pom[2]);
        break;
    case 26:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &tags[5], &pom[0], &pom[1], &pom[2]);
        break;
    case 27:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &tags[5], &tags[6], &pom[0], &pom[1], &pom[2]);
        break;
    case 31:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 32:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 33:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 34:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 35:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 36:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &tags[5], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 37:
        nred=sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &tags[2], &tags[3], &tags[4], &tags[5], &tags[6], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    
    // 3D Elements
    case 42: // Tetrahedron (4 nodes)
        nred = sscanf(buf, "%d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2], &pom[3]);
        break;
    case 52: // Hexahedron (8 nodes)
        nred = sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2], &pom[3], &pom[4], &pom[5], &pom[6], &pom[7]);
        break;
    case 62: // Prism/Wedge (6 nodes)
        nred = sscanf(buf, "%d %d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2], &pom[3], &pom[4], &pom[5]);
        break;
    case 72: // Pyramid (5 nodes)
        nred = sscanf(buf, "%d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2], &pom[3], &pom[4]);
        break;

    default:
        writelog("Error[MSH]:\t Incorrect number of tags (NTAGS %d)\n", tgs);
        break;
    }
    p->tagPhysical    = tags[0];
    p->tagGeometrical = tags[1];
    p->tagSubDomain   = 1;
    if (tgs > 3)
	p->tagSubDomain   = tags[3];

    // printf("tgs: %d \t nred: %d \t val: %d \t etp: %d \t idx: %d\n",tgs, nred, val, etp, idx);

    switch (etp)
    {
    case GMSH_SEGMENT:
        if (nred != tgs + 3 + 2)
            writelog("Error[MSH]:\t Incorrect mesh format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        break;
    case GMSH_TRIANGLE:
        if (nred != tgs + 3 + 3)
            writelog("Error[MSH]:\t Incorrect mesh format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        p->ilist[2] = pom[2] - 1;
        break;
    case GMSH_QUADRANGLE:
        if (nred != tgs + 3 + 4)
            writelog("Error[MSH]:\t Incorrect mesh format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        p->ilist[2] = pom[2] - 1;
        p->ilist[3] = pom[3] - 1;
        break;
    case GMSH_TETRAHEDRON: // Tetrahedron element with 4 nodes (3D)
        if (nred != tgs + 3 + 4) 
            writelog("Error[MSH]: Incorrect TETRAHEDRON format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        p->ilist[2] = pom[2] - 1;
        p->ilist[3] = pom[3] - 1;
        break;
    case GMSH_HEXAHEDRON: // Hexahedron element with 8 nodes (3D)
        if (nred != tgs + 3 + 8) 
            writelog("Error[MSH]: Incorrect HEXAHEDRON format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        p->ilist[2] = pom[2] - 1;
        p->ilist[3] = pom[3] - 1;
        p->ilist[4] = pom[4] - 1;
        p->ilist[5] = pom[5] - 1;
        p->ilist[6] = pom[6] - 1;
        p->ilist[7] = pom[7] - 1;
        break;
    case GMSH_WEDGE: // Wedge (prism) element with 6 nodes (3D)
        if (nred != tgs + 3 + 6) 
            writelog("Error[MSH]: Incorrect WEDGE format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        p->ilist[2] = pom[2] - 1;
        p->ilist[3] = pom[3] - 1;
        p->ilist[4] = pom[4] - 1;
        p->ilist[5] = pom[5] - 1;
        break;
    case GMSH_PYRAMID: // Pyramid element with 5 nodes (3D)
        // nred = sscanf(buf, "%d %d %d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2], &pom[3], &pom[4]);
        if (nred != tgs + 3 + 5)
            writelog("Error[MSH]: Incorrect PYRAMID format.\n");
        p->ilist[0] = pom[0] - 1;
        p->ilist[1] = pom[1] - 1;
        p->ilist[2] = pom[2] - 1;
        p->ilist[3] = pom[3] - 1;
        p->ilist[4] = pom[4] - 1;
        break;
    }

    return 0;
}
/* ------------------------------------------------------------------------------------------------------- */
void   gmshfile_add2physical(gmshfile *p, int mark, int etyp)
{
    int i; 

    if (mark == 0)
        return;
    for(i = 0; i < p->nPhysical; i++)
        if (p->labelPhysical[i] == mark)
            break;
    if (i == p->nPhysical)
    {
        if (i >= MAX_PHYSICAL)
            CriticalErrorMsg("Error[GMSH]:\t Number of subdomains too large!\n");
        p->labelPhysical[i] = mark;
        p->ntriPhysical[i]  = 0;
        p->nquadPhysical[i] = 0;
        p->nedgePhysical[i] = 0;
        p->ntetPhysical[i]  = 0; // Added for tetrahedra
        p->nhexPhysical[i]  = 0; // Added for hexahedra
        p->nwedgePhysical[i] = 0; // Added for wedges
        p->npyrPhysical[i] = 0; // Added for pyramid
        p->nPhysical++;
    }
  
    switch (etyp)
    {
    case GMSH_SEGMENT:     p->nedgePhysical[i]++; break;
    case GMSH_TRIANGLE:    p->ntriPhysical[i]++;  break;
    case GMSH_QUADRANGLE:  p->nquadPhysical[i]++; break;
    case GMSH_TETRAHEDRON:  p->ntetPhysical[i]++;      break; // 3D
    case GMSH_HEXAHEDRON:   p->nhexPhysical[i]++;      break; // 3D
    case GMSH_WEDGE:        p->nwedgePhysical[i]++;    break; // 3D
    case GMSH_PYRAMID:      p->npyrPhysical[i]++;      break; // 3D
    default:
        CriticalErrorMsg("GMSH add physical domain attributed failed, etyp %d\n", etyp);
        break;
    }
}
/* ------------------------------------------------------------------------------------------------------- */
void   gmshfile_add2subdomain(gmshfile *p, int id, int etyp)
{
  int i; 
  if (id == 0)
    return;
  for(i = 0; i < p->nSubdomain; i++)
    if (p->labelSubdomain[i] == id)
      break;
  
  if (i == p->nSubdomain)
  {
      if (i >= MAX_SUBDOMAINS)
          CriticalErrorMsg("Error[GMSH]:\t Number of subdomains too large!\n");
      p->labelSubdomain[i] = id;
      p->ntriSubdomain[i]  = 0;
      p->nquadSubdomain[i] = 0;
      p->nedgeSubdomain[i] = 0;
      p->ntetSubdomain[i]  = 0;  // Added for tetrahedra
      p->nhexSubdomain[i]  = 0;  // Added for hexahedra
      p->nwedgeSubdomain[i] = 0; // Added for wedges
      p->npyrSubdomain[i] = 0; // Added for pyramid
      p->nSubdomain++;
  }
  switch (etyp)
  {
  case GMSH_SEGMENT:      p->nedgeSubdomain[i]++; break;
  case GMSH_TRIANGLE:     p->ntriSubdomain[i]++;  break;
  case GMSH_QUADRANGLE:   p->nquadSubdomain[i]++; break;
  case GMSH_TETRAHEDRON:  p->ntetSubdomain[i]++;      break; // 3D
  case GMSH_HEXAHEDRON:   p->nhexSubdomain[i]++;      break; // 3D
  case GMSH_WEDGE:        p->nwedgeSubdomain[i]++;    break; // 3D
  case GMSH_PYRAMID:      p->npyrSubdomain[i]++;      break; // 3D
  default: 
      CriticalErrorMsg("GMSH add subdomain attributed failed, etyp %d\n", etyp);
      break;
  }
  
}
/* ------------------------------------------------------------------------------------------------------- */
void   gmshfile_report(gmshfile *p, const char *fname)
{
    int ii;
    writelog("Init[GMSH]: \t\t GMSH file %s statistics: NP %d ", fname, p->NPoints);
    writelog("NE %d NT %d NQ %d ", p->NElements, p->NTriangles, p->NQuads);
    writelog("NS %d NBS %d\n", p->NEdges, p->NBndrEdges);

    writelog("           \t\t NPHYSICAL   %d", p->nPhysical);
    for(ii = 0; ii < p->nPhysical; ii++)
	if (p->nPhysical > 6)
	    writelog("\n           \t\t (MARK %d, NT %d NQ %d NS %d)", p->labelPhysical[ii], p->ntriPhysical[ii], p->nquadPhysical[ii], p->nedgePhysical[ii]);
	else
	{
	    writelog(" (M %d, NT %d NQ %d NS %d)", p->labelPhysical[ii], p->ntriPhysical[ii], p->nquadPhysical[ii], p->nedgePhysical[ii]);
	}
    writelog("\n           \t\t NSUBDOMAINS %d", p->nSubdomain);
    for(ii = 0; ii < p->nSubdomain; ii++)
	if (p->nSubdomain > 6)
	    writelog("           \t\t (MARK %d, NT %d NQ %d NS %d)\n", p->labelSubdomain[ii], p->ntriSubdomain[ii], p->nquadSubdomain[ii], p->nedgeSubdomain[ii]);
	else
	    writelog(" (M %d, NT %d NQ %d NS %d)", p->labelSubdomain[ii], p->ntriSubdomain[ii], p->nquadSubdomain[ii], p->nedgeSubdomain[ii]);
	    
    writelog("\n");
}
/* ------------------------------------------------------------------------------------------------------- */
void   gmshfile_read(gmshfile *p, const char *fname)
{
    int      i;
    FILE     *fid;
    char     buf[BUFFER_SIZE];
    gmshline current;
    int dim = 2;
    
    p->NPoints       = 0;
    p->NGeomElements = 0;
    p->NElements     = 0;
    p->NTriangles    = 0;
    p->NQuads        = 0;
    p->NTetrahedra   = 0; // New for 3D
    p->NHexahedra    = 0; // New for 3D
    p->NWedges       = 0; // New for 3D
    p->NPyramids     = 0; // New for 3D
    p->nSubdomain    = 0;
    p->nPhysical     = 0;
    p->NBndrEdges    = 0;

    current.tagPhysical = 0;
    current.etyp     = 0;

    fid = fopen(fname, "r"); FileCheck(fid, "gmshfile_read opening %s\n", fname);
    while (!feof(fid))
    { 
        fgets(buf, BUFFER_SIZE, fid);
        if (memcmp(buf, "$Nodes", 6) == 0)
        {  
            fgets(buf, BUFFER_SIZE, fid);
            sscanf(buf, "%d", &p->NPoints);
        }

        if (memcmp(buf, "$Elements", 9) == 0)
        {
            fgets(buf, BUFFER_SIZE, fid);
            sscanf(buf, "%d", &p->NGeomElements);
            for(i = 0; i < p->NGeomElements; i++)
            {
                fgets(buf, BUFFER_SIZE, fid);
                gmshline_read(&current, buf);
                gmshfile_add2subdomain(p, current.tagSubDomain ,current.etyp);
                gmshfile_add2physical( p, current.tagPhysical  ,current.etyp);
                switch (current.etyp)
                {
                case GMSH_SEGMENT:
                    if (current.tagPhysical < 0) 
                        p->NInterfaceEdges++;
                    else
                        p->NBndrEdges++; 
                    break;
                case GMSH_TRIANGLE:
                    p->NTriangles++; 
                    dim = 2;
                    break;
                case GMSH_QUADRANGLE:
                    p->NQuads++; 
                    dim = 2;
                    break;

                case GMSH_TETRAHEDRON: // New for 3D
                    p->NTetrahedra++;
                    dim = 3;
                    break;
                case GMSH_HEXAHEDRON: // New for 3D
                    p->NHexahedra++;
                    dim = 3;
                    break;
                case GMSH_WEDGE: // New for 3D
                    p->NWedges++;
                    dim = 3;
                    break;
                case GMSH_PYRAMID: // New for 3D
                    p->NPyramids++;
                    dim = 3;
                    break;

                default:
                    writelog("\n\nError[MSH]:\t Unknown element %d, type %d!\n", i + 1, current.etyp);
                    exit(1);
                    break;
                }
            }
        }
    }
    fclose(fid);

    // p->NEdges         = (p->NTriangles * 3 + p->NQuads * 4 + p->NBndrEdges) / 2 ;
    // if ( (p->NTriangles * 3 + p->NQuads * 4 + p->NBndrEdges) % 2 )
    //     CriticalErrorMsg("Error[MSH]:\t Something is wrong, mesh is not consistent! Remove internal boundaries! \n");
    
    if (dim == 2)
    {
        p->NEdges = (p->NTriangles * 3 + p->NQuads * 4 + p->NBndrEdges) / 2;
        if ((p->NTriangles * 3 + p->NQuads * 4 + p->NBndrEdges) % 2)
            CriticalErrorMsg("Error[MSH]: Mesh is inconsistent in 2D.\n");
    }
    else if (dim == 3)
    {
       p->NEdges = 0;
    }

    gmshfile_report(p, fname);
}
/* ------------------------------------------------------------------------------------------------------- */
void   gmshfile_read4multiphysics(gmshfile *p, const char *fname, short isInteriorEdgeMark(short))
{
    int      i;
    FILE     *fid;
    char     buf[BUFFER_SIZE];
    gmshline current;
    
    p->NPoints       = 0;
    p->NGeomElements = 0;
    p->NElements     = 0;
    p->NTriangles    = 0;
    p->NQuads        = 0;
    p->nSubdomain    = 0;
    p->nPhysical     = 0;
    p->NBndrEdges    = 0;
    current.tagPhysical = 0;
    current.etyp     = 0;

    fid = fopen(fname, "r"); FileCheck(fid, "gmshfile_read opening %s\n", fname);
    while (!feof(fid)) 
    { 
        fgets(buf, BUFFER_SIZE, fid); 
        if (memcmp(buf, "$Nodes", 6) == 0)
        {  
            fgets(buf, BUFFER_SIZE, fid);
            sscanf(buf, "%d", &p->NPoints);
        }

        if (memcmp(buf, "$Elements", 9) == 0)
        {
            fgets(buf, BUFFER_SIZE, fid);
            sscanf(buf, "%d", &p->NGeomElements);
            for(i = 0; i < p->NGeomElements; i++)
            {
                fgets(buf, BUFFER_SIZE, fid);
                gmshline_read(&current, buf);
                gmshfile_add2subdomain(p, current.tagSubDomain ,current.etyp);
                gmshfile_add2physical( p, current.tagPhysical  ,current.etyp);
                switch (current.etyp)
                {
                case GMSH_SEGMENT:
                    if (isInteriorEdgeMark(current.tagPhysical)) 
                        p->NInterfaceEdges++;
                    else
                        p->NBndrEdges++; 
                    break;
                case GMSH_TRIANGLE:
                    p->NTriangles++; 
                    break;
                case GMSH_QUADRANGLE:
                    p->NQuads++; 
                    break;
                default:
                    writelog("\n\nError[MSH]:\t Unknown element %d, type %d!\n", i + 1, current.etyp);
                    exit(1);
                    break;
                }
            }
        }
    }
    fclose(fid);
    p->NEdges         = (p->NTriangles * 3 + p->NQuads * 4 + p->NBndrEdges) / 2 ;
    if ( (p->NTriangles * 3 + p->NQuads * 4 + p->NBndrEdges) % 2 )
        CriticalErrorMsg("Error[MSH]:\t Something is wrong, mesh is not consistent!\n");
    gmshfile_report(p, fname);
}
/* ------------------------------------------------------------------------------------------ */
int    gmsh_getNPoints(const char *fname)
{
    FILE *fid;
    char buf[BUFFER_SIZE];
    int  numPoints;
    
    debuglog("Debug[GMSH]:\t\t gmsh::getNPoints\n");
    numPoints = 0;
    fid = fopen(fname, "r"); FileCheck(fid, "gmsh_getNPoints, file %s\n", fname);
    while (!feof(fid)) 
    { 
        fgets(buf, BUFFER_SIZE, fid); 
        if( memcmp(buf, "$Nodes", 6) == 0)
        {  
            fgets(buf, BUFFER_SIZE, fid);
            sscanf(buf, "%d", &numPoints);
            break;
        }
    }
    fclose(fid);
    debuglog("Debug[GMSH]:\t\t Done.\n");
    return numPoints;
}
/* ------------------------------------------------------------------------------------------ */
int    gmsh_getNodesCoordinates(FILE *fid, double *x, double *y, double *z, int dim)
{
    char buf[BUFFER_SIZE];
    int i, np, idx, ok;
    double _xx, _yy, _zz;

    debuglog("Debug[GMSH]:\t\t Load xy-coordinates of vertices from file.\n");
    do 
        fgets(buf, BUFFER_SIZE, fid); 
    while ((strncasecmp(buf, "$Nodes", 6) != 0) && (feof(fid) == 0) );
    // printf("dim: %d\n", dim);
    fgets(buf, BUFFER_SIZE, fid); 
    sscanf(buf, "%d", &np);
    for(i = 0; i < np; i++)
    {
        fgets(buf, BUFFER_SIZE, fid);

        if (dim == 2)
        {
            ok = sscanf(buf, "%d %lf %lf", &idx, &_xx, &_yy);
            x[idx - 1] = _xx;
            y[idx - 1] = _yy;
        }
        else if (dim == 3) 
        {
            ok = sscanf(buf, "%d %lf %lf %lf", &idx, &_xx, &_yy, &_zz);
            x[idx - 1] = _xx;
            y[idx - 1] = _yy;
            z[idx - 1] = _zz;
            // printf("idx: %d x: %g y: %g z: %g\n", idx-1, _xx, _yy, _zz);
        }
        
        if ( (dim == 2 && ok != 3) || (dim == 3 && ok != 4) )
            CriticalErrorMsg("gmsh::getNodesCoordinates failed, probably incorrect mesh format\n");
    }

    fgets(buf, BUFFER_SIZE, fid);
    if (strncasecmp(buf, "$EndNodes", 9) != 0)
        CriticalErrorMsg("gmsh::getNodesCoordinates failed, probably incorrect mesh format (EndNodes - missing)\n");
    debuglog("Debug[GMSH]:\t\t Load xy-coordinates of vertices \t\t(done)\n");
    return np;

}
