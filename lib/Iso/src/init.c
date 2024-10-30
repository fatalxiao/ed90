#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(pava)(double *y, double *w, int *kt, int *n);
extern void F77_NAME(smooth)(int *nrow, int *ncol, int *ndim, double *x,
                     double *w, double *a, double *b, int *ncycle, int *icycle,
                     double *g, double *eps1, double *eps2, int *ifault, double *fx,
                     double *pw, double *w1, double *wt, int *nw);
extern void F77_NAME(ufit)(double *xk, double *wk, double *xmode, double *x,
	       	double *w, double *mse, double *x1, double *w1, double *x2,
	       	double *w2, int *ind, int *kt, int *n);
extern void F77_NAME(unimode)(double *y, double *w, double *y1, double *w1,
	       	double *y2, double *w2, int *ind, int *kt, double *tau, int *n);

/*
   Note that the unimode routine does not feature in the foregoing
   since unimode is called only by ufit and never called directly
   by .Fortran().
*/

static const R_FortranMethodDef FortranEntries[] = {
    {"pava",   (DL_FUNC) &F77_NAME(pava),    4},
    {"smooth", (DL_FUNC) &F77_NAME(smooth), 18},
    {"ufit",   (DL_FUNC) &F77_NAME(ufit),   13},
    {"unimode", (DL_FUNC) &F77_NAME(unimode),10},
    {NULL, NULL, 0}
};

void R_init_Iso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
