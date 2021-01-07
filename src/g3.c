#include <math.h>
#include <R.h>
#include "geom3.h"
#include "functable.h"

/*
	$Revision: 1.3 $	$Date: 2012/05/22 07:17:31 $

	G function (nearest neighbour distribution) of 3D point pattern


	Let 	b = distance from point p[i] to boundary of box
	 	d = distance from p[i] to nearest p[j] 


	method = 1	naive ratio estimator (Ripley 1981)

			numerator(r)  = count(i: b >= r, d <= r)
			denominator(r)  = count(i: b >= r)

	method = 2	minus sampling estimator

			numerator(r) = count(i: b >= r, d <= r)
			denominator(r) = lambda * volume(x: b >= r)

			where lambda = (no of points)/volume(box)

	method = 3	Hanisch's G3

			numerator(r) = count(i: b >= d, d <= r)
			denominator(r) = count(i: b >= d)

	method = 4	Hanisch's G4

			numerator(r) = count(i: b >= d, d <= r)
			denominator(r) = fudge * volume(x: b >= r)

			fudge = numerator(R)/denominator(R)
			R = sup{r : denominator(r) > 0 }


# /////////////////////////////////////////////
# AUTHOR: Adrian Baddeley, CWI, Amsterdam, 1991.
#
# MODIFIED BY: Adrian Baddeley, Perth 2009, 2012.
#
# This software is distributed free
# under the conditions that
# 	(1) it shall not be incorporated
# 	in software that is subsequently sold
# 	(2) the authorship of the software shall
# 	be acknowledged in any publication that 
# 	uses results generated by the software
# 	(3) this notice shall remain in place
# 	in each file.
# //////////////////////////////////////////////


*/

#define MIN(X,Y) (((X) > (Y)) ? (Y) : (X))

double *
nndist3(p, n, b)
		/* compute nearest neighbour distance for each p[i] */
     Point *p;
     int n;
     Box *b;
{
  register int i, j;
  register double dx, dy, dz, dist2, nearest2, huge2;
  Point *ip, *jp;
  double *nnd;

  nnd = (double *) R_alloc(n, sizeof(double));

  dx = b->x1 - b->x0;
  dy = b->y1 - b->y0;
  dz = b->z1 - b->z0;
  huge2 = 2.0 * (dx * dx + dy * dy + dz * dz);
	
  /* scan each point and find closest */
  for( i = 0; i < n; i++) {
    ip = p + i;
    nearest2 = huge2;
    for(j = 0; j < n; j++)
      if(j != i) {
	jp = p + j;
	dx = jp->x - ip->x;
	dy = jp->y - ip->y;
	dz = jp->z - ip->z;
	dist2 = dx * dx + dy * dy + dz * dz;
	if(dist2 < nearest2)
	  nearest2 = dist2;
      }
    nnd[i] = sqrt(nearest2);
  }
  return(nnd);
}

double *
border3(p, n, b)
		/* compute distances to border */
     Point *p;
     int n;
     Box *b;
{
  register int i;
  register double bord;
  register Point *ip;
  double *bored;

  bored = (double *) R_alloc(n, sizeof(double));

  for( i = 0; i < n; i++) {
    ip = p + i;
    bord = MIN(ip->x - b->x0, b->x1 - ip->x);
    bord = MIN(bord, ip->y - b->y0);
    bord = MIN(bord, b->y1 - ip->y);
    bord = MIN(bord, ip->z - b->z0);
    bord = MIN(bord, b->z1 - ip->z);
    bored[i] = bord;
  }
  return(bored);
}

void
g3one(p, n, b, g)
     Point *p;
     int n;
     Box *b;
     Ftable *g;
{
  register int i, l, lbord, lnnd;
  double dt;
  double	*bord, *nnd;

  bord = border3(p, n, b);
  nnd  = nndist3(p, n, b);
	
  /* initialise */
  for(l = 0; l < g->n; l++)
    (g->num)[l] = (g->denom)[l] = 0.0;

  /* spacing of argument in result vector g */
  dt = (g->t1 - g->t0)/(g->n - 1);

  for(i = 0; i < n; i++) {
    lbord = floor( (bord[i] - g->t0) / dt );
    if(lbord >= g->n) 
      lbord = g->n - 1;
    for(l = 0; l <= lbord; l++)
      (g->denom)[l] += 1.0;

    lnnd = ceil( (nnd[i] - g->t0) / dt );
    if(lnnd < 0) lnnd = 0;
    for(l = lnnd; l <= lbord; l++)
      (g->num)[l] += 1.0;
  }

  /* compute ratio */
  for(l = 0; l < g->n; l++)
    (g->f)[l] = ((g->denom)[l] > 0)?
      (g->num)[l] / (g->denom)[l] : 1.0;
			   
}

void
g3three(p, n, b, g)
     Point *p;
     int n;
     Box *b;
     Ftable *g;
{
  register int i, l, lmin;
  double dt;
  int	denom;
  double	*bord, *nnd;

  bord = border3(p, n, b);
  nnd  = nndist3(p, n, b);
	
  /* initialise */
  denom = 0;
  for(l = 0; l < g->n; l++)
    (g->num)[l]   = 0.0;
  
  /* spacing of argument in result vector g */
  dt = (g->t1 - g->t0)/(g->n - 1);

  for(i = 0; i < n; i++) { 
    if(nnd[i] <= bord[i]) {
      ++denom;

      lmin = ceil( (nnd[i] - g->t0) / dt );
      if(lmin < 0) lmin = 0;
      for(l = lmin; l < g->n; l++)
	(g->num)[l] += 1.0;
    }
  }

  /* compute ratio */
  for(l = 0; l < g->n; l++) {
    (g->denom)[l] = denom;
    (g->f)[l] = (denom > 0)?
      (g->num)[l] / (double) denom
      : 1.0;
  }
}

void
g3cen(p, n, b, count)
     Point *p;
     int n;
     Box *b;
     H4table *count;
{
  register int i, lcen, lobs;
  register double dt, cens, obsv;
  double	*bord, *nnd;

  bord = border3(p, n, b);
  nnd  = nndist3(p, n, b);
	
  /* spacing of histogram cells */
  dt = (count->t1 - count->t0)/(count->n - 1);

  /* 'count' is assumed to have been initialised */
  for(i = 0; i < n; i++) { 
    obsv = nnd[i];
    cens = bord[i];
    lobs = ceil( (obsv - count->t0) / dt );
    lcen = floor( (cens - count->t0) / dt );
    if(obsv <= cens) {
      /* observation is uncensored; 
	 increment all four histograms */
      if(lobs >= count->n)
	++(count->upperobs);
      else if(lobs >= 0) {
	(count->obs)[lobs]++;
	(count->nco)[lobs]++;
      }
      if(lcen >= count->n)
	++(count->uppercen);
      else if(lcen >= 0) {
	(count->cen)[lcen]++;
	(count->ncc)[lcen]++;
      }
    } else {
      /* observation is censored; 
	 increment only two histograms */
      lobs = MIN(lobs, lcen);
      if(lobs >= count->n)
	++(count->upperobs);
      else if(lobs >= 0) 
	(count->obs)[lobs]++;

      if(lcen >= count->n)
	++(count->uppercen);
      else if(lcen >= 0) 
	(count->cen)[lcen]++;
    }
  }
}
