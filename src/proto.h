#include <R.h>
#include <Rinternals.h>

/*
  Prototype declarations for all native routines in spatstat.core package

  Automatically generated - do not edit! 

*/

/*
  
                  Functions invoked by .C

*/

void delta2area(double *, double *, double *, double *, int *, double *, double *, double *, double *, int *); 
void delta2area(double *, double *, double *, double *, int *, double *, double *, double *, double *, int *);
void digberJ(double *, double *, int *, int *, int *, double *);
void Gdenspt(int *, double *, double *, double *, double *); 
void Gwtdenspt(int *, double *, double *, double *, double *, double *); 
void Gwtdenspt(int *, double *, double *, double *, double *, double *); 
void denspt(int *, double *, double *, double *, double *, double *); 
void wtdenspt(int *, double *, double *, double *, double *, double *, double *); 
void wtdenspt(int *, double *, double *, double *, double *, double *, double *); 
void adenspt(int *, double *, double *, double *, double *, double *, double *); 
void awtdenspt(int *, double *, double *, double *, double *, double *, double *, double *); 
void awtdenspt(int *, double *, double *, double *, double *, double *, double *, double *); 
void crdenspt(int *, double *, double *, int *, double *, double *, double *, double *, double *); 
void wtcrdenspt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void wtcrdenspt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void acrdenspt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void awtcrdenspt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *); 
void awtcrdenspt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *);
void segdens(double *, int *, double *, double *, double *, double *, int *, double *, double *, double *); 
void segwdens(double *, int *, double *, double *, double *, double *, double *, int *, double *, double *, double *);
void Ediggra(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *, double *); 
void ESdiggra(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *, double *, int *);
void Ediggatsti(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *);
void ripleybox(int *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void ripleypoly(int *, double *, double *, double *, int *, double *, int *, double *, double *, double *, double *, double *); 
void rippolDebug(int *, double *, double *, double *, int *, double *, int *, double *, double *, double *, double *, double *);
void RcallK3(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, int *, double *, double *, double *, int *); 
void RcallG3(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, int *, double *, double *, double *, int *); 
void RcallF3(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *); 
void RcallF3cen(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *); 
void RcallG3cen(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *); 
void Rcallpcf3(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, int *, double *, double *, double *, int *, double *); 
void RcallF3(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *);
void locxprod(int *, double *, double *, int *, double *, double *, double *, int *, double *, double *);
void Efiksel(int *, double *, double *, int *, double *, double *, double *, double *, double *);
void Egeyer(int *, double *, double *, int *, int *, double *, double *, int *, double *, double *, double *);
void Cidw(double *, double *, double *, int *, double *, double *, int *, double *, double *, int *, double *, double *, double *, double *); 
void Cidw2(double *, double *, double *, int *, double *, double *, int *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void idwloo(double *, double *, double *, int *, double *, double *, double *, double *); 
void idwloo2(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *);
void locprod(int *, double *, double *, double *, int *, double *, double *); 
void locxprod(int *, double *, double *, int *, double *, double *, double *, int *, double *, double *);
void KborderI(int *, double *, double *, double *, int *, double *, int *, int *); 
void KborderD(int *, double *, double *, double *, int *, double *, double *, double *); 
void Kwborder(int *, double *, double *, double *, double *, int *, double *, double *, double *); 
void KnoneI(int *, double *, double *, int *, double *, int *); 
void KnoneD(int *, double *, double *, int *, double *, double *); 
void Kwnone(int *, double *, double *, double *, int *, double *, double *); 
void KrectWtd(double *, double *, int *, double *, double *, double *, int *, double *, double *, int *, int *, int *, int *, double *, double *, double *, double *, double *); 
void KrectInt(double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, int *, int *, double *, double *, int *, int *, int *); 
void KrectDbl(double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, int *, int *, double *, double *, double *, double *, double *);
void locpcfx(int *, double *, double *, int *, int *, double *, double *, int *, int *, double *, double *, double *); 
void locWpcfx(int *, double *, double *, int *, int *, double *, double *, int *, double *, int *, double *, double *, double *);
void knownCif(char *, int *);
void scantrans(double *, double *, int *, double *, double *, double *, double *, int *, int *, double *, int *);
void Gsmoopt(int *, double *, double *, double *, int *, double *, double *); 
void Gwtsmoopt(int *, double *, double *, double *, int *, double *, double *, double *); 
void smoopt(int *, double *, double *, double *, int *, double *, double *, double *); 
void wtsmoopt(int *, double *, double *, double *, int *, double *, double *, double *, double *); 
void asmoopt(int *, double *, double *, double *, int *, double *, double *, double *); 
void awtsmoopt(int *, double *, double *, double *, int *, double *, double *, double *, double *); 
void crsmoopt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void wtcrsmoopt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *); 
void acrsmoopt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *); 
void awtcrsmoopt(int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *);
void Cclosepaircounts(int *, double *, double *, double *, int *); 
void Ccrosspaircounts(int *, double *, double *, int *, double *, double *, double *, int *);
/*

             Functions invoked by .Call

*/
SEXP thinjumpequal(SEXP, SEXP, SEXP);
SEXP xmethas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP xmethas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP PerfectStrauss(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP PerfectHardcore(SEXP, SEXP, SEXP, SEXP); 
SEXP PerfectStraussHard(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP PerfectDiggleGratton(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP PerfectDGS(SEXP, SEXP, SEXP, SEXP); 
SEXP PerfectPenttinen(SEXP, SEXP, SEXP, SEXP, SEXP);
