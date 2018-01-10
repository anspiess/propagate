#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/
  
/* .Call calls */
extern SEXP C_colVarsC(SEXP);
extern SEXP C_rowVarsC(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_colVarsC", (DL_FUNC) &C_colVarsC, 1},
  {"C_rowVarsC", (DL_FUNC) &C_rowVarsC, 1},
  {NULL, NULL, 0}
};

void R_init_propagate(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
