#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
/*extern void davies(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);*/
extern void davies(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lq, int *lim1, double *acc, double* trace, int* ifault, double *res);

static const R_CMethodDef CEntries[] = {
    {"davies", (DL_FUNC) &davies, 13},
    {NULL, NULL, 0}
};

void R_init_DiceOptim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
