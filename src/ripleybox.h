/*

  ripleybox.h

  Ripley's edge correction for rectangular windows

  This file is #included multiple times in corrections.c
  Macros used:

  RIPLEYFUN      Name of C function
  DEBUGBOX       #defined if debugging information should be printed.

  *CHUNKLOOP     defined in chunkloop.h
  TWOPI          defined in Rmath.h

  $Revision: 1.3 $     $Date: 2021/10/31 06:40:58 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
  Licence: GNU Public Licence >= 2

*/

void RIPLEYFUN(nx, x, y, rmat, nr, xmin, ymin, xmax, ymax,  epsilon, out)
     /* inputs */
     int *nx, *nr;  /* dimensions */
     double *x, *y; /* coordinate vectors of length nx */
     double *rmat;  /* matrix nx by nr  */
     double *xmin, *ymin, *xmax, *ymax;  /* box dimensions */
     double *epsilon; /* threshold for proximity to corner */
     /* output */
     double *out;  /* output matrix nx by nr */
{
  int i, j, n, m, ijpos, ncor, maxchunk;
  double xx, yy, x0, y0, x1, y1, dL, dR, dU, dD, aL, aU, aD, aR, rij;
  double cL, cU, cD, cR, bLU, bLD, bRU, bRD, bUL, bUR, bDL, bDR;
  double corner, extang;
  double eps;

  n  = *nx;
  m  = *nr;
  x0 = *xmin;
  y0 = *ymin;
  x1 = *xmax;
  y1 = *ymax;
  eps = *epsilon;

  OUTERCHUNKLOOP(i, n, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 16384) {
      xx = x[i];
      yy = y[i];
      /* 
	 perpendicular distance from point to each edge of rectangle
	 L = left, R = right, D = down, U = up
      */
      dL = xx - x0;
      dR = x1 - xx;
      dD = yy - y0;
      dU = y1 - yy;

      /*
	test for corner of the rectangle
      */
#define ABS(X) (((X) >= 0) ? (X) : (-X))
#define SMALL(X) ((ABS(X) < eps) ? 1 : 0)

      ncor = SMALL(dL) + SMALL(dR) + SMALL(dD) + SMALL(dU);
      corner = (ncor >= 2) ? YES : NO;
  
      /* 
	 angle between 
	 - perpendicular to edge of rectangle
	 and 
	 - line from point to corner of rectangle

      */
      bLU = atan2(dU, dL);
      bLD = atan2(dD, dL);
      bRU = atan2(dU, dR);
      bRD = atan2(dD, dR);
      bUL = atan2(dL, dU);
      bUR = atan2(dR, dU);
      bDL = atan2(dL, dD);
      bDR = atan2(dR, dD);

      for(j = 0; j < m; j++) {
	ijpos = j * n + i;
	rij = rmat[ijpos];
#ifdef DEBUGBOX
	Rprintf("rij = %lf\n", rij);
#endif
	if(rij == 0.0) {
	  /* Circle of radius 0 */
	  out[ijpos] = 1.0;
	} else {
	  /*
	    Fraction of circle
	    Compute half the angle subtended by the intersection between
	    the circle of radius r[i,j] centred on point i
	    and each edge of the rectangle (prolonged to an infinite line)
	  */
	  aL = (dL < rij) ? acos(dL/rij) : 0.0;
	  aR = (dR < rij) ? acos(dR/rij) : 0.0;
	  aD = (dD < rij) ? acos(dD/rij) : 0.0;
	  aU = (dU < rij) ? acos(dU/rij) : 0.0;
#ifdef DEBUGBOX
	  Rprintf("aL = %lf\n", aL);
	  Rprintf("aR = %lf\n", aR);
	  Rprintf("aD = %lf\n", aD);
	  Rprintf("aU = %lf\n", aU);
#endif
	  /* apply maxima */

	  cL = MIN(aL, bLU) + MIN(aL, bLD);
	  cR = MIN(aR, bRU) + MIN(aR, bRD);
	  cU = MIN(aU, bUL) + MIN(aU, bUR);
	  cD = MIN(aD, bDL) + MIN(aD, bDR);
#ifdef DEBUGBOX
	  Rprintf("cL = %lf\n", cL);
	  Rprintf("cR = %lf\n", cR);
	  Rprintf("cD = %lf\n", cD);
	  Rprintf("cU = %lf\n", cU);
#endif

	  /* total exterior angle over 2 pi */
	  extang = (cL + cR + cU + cD)/TWOPI;

#ifdef DEBUGBOX
	  Rprintf("extang = %lf\n", extang);
#endif

	  /* add pi/2 for corners */
	  if(corner) {
	    extang += 1.0/4.0;
#ifdef DEBUGBOX
	    Rprintf("extang = %lf\n", extang);
#endif
	  }
	
	  /* OK, now compute weight */
	  out[ijpos] = 1 / (1 - extang);
	}
      }
    }
  }
}

