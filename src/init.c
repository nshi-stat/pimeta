#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _pimeta_pwchisqCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pimeta_bootPICppWrap(SEXP, SEXP, SEXP, SEXP);
extern SEXP _pimeta_rtau2CppWrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_pimeta_pwchisqCpp", (DL_FUNC) &_pimeta_pwchisqCpp,  8},
  {"_pimeta_bootPICppWrap", (DL_FUNC) &_pimeta_bootPICppWrap,  4},
  {"_pimeta_rtau2CppWrap",  (DL_FUNC) &_pimeta_rtau2CppWrap,  11},
  {NULL, NULL, 0}
};

void R_init_pimeta(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}