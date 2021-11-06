/*

  corrections.c

  Edge corrections

  $Revision: 1.17 $     $Date: 2021/10/25 10:18:31 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"
#include "yesno.h"
#include "constants.h"

/* This constant is defined in Rmath.h */
#define TWOPI M_2PI

#define MIN(A,B) (((A) < (B)) ? (A) : (B))

#define BETWEEN(X,X0,X1) ((double) ( ( (X) - (X0) ) * ( (X) - (X1) ) ) <= 0.0)

#define UNDER(X,Y,X0,Y0,X1,Y1) \
  ((double) ( ( (Y1) - (Y0) ) * ( (X) - (X0) ) ) >= (double) ( ( (Y) - (Y0) ) * ( (X1) - (X0) ) ) )

#define UNDERNEATH(X,Y,X0,Y0,X1,Y1) \
  ((((double) (X0)) < ((double) (X1))) ? UNDER(X,Y,X0,Y0,X1,Y1) : UNDER(X,Y,X1,Y1,X0,Y0))

#define TESTINSIDE(X,Y,X0,Y0,X1,Y1) \
  (BETWEEN(X,X0,X1) && UNDERNEATH(X, Y, X0, Y0, X1, Y1))



/* C function ripleybox */
#undef DEBUGBOX
#define RIPLEYFUN ripleybox
#include "ripleybox.h"
#undef RIPLEYFUN

/* C function ripboxDebug */
#define DEBUGBOX
#define RIPLEYFUN ripboxDebug
#include "ripleybox.h"
#undef RIPLEYFUN
#undef DEBUGBOX

/* C function ripleypoly */
#undef DEBUGPOLY
#define RIPLEYFUN ripleypoly
#include "ripleypoly.h"
#undef RIPLEYFUN

/* C function rippolDebug */
#define DEBUGPOLY
#define RIPLEYFUN rippolDebug
#include "ripleypoly.h"
#undef RIPLEYFUN
#undef DEBUGPOLY


