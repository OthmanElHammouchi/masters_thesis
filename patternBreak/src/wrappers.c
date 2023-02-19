#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rinterface.h>

void F77_NAME(mack_boot)(int n_dev, double triangle[n_dev][n_dev], int resids_type, int boot_type, int dist, int n_boot, double reserve[n_boot]);

SEXP mack_boot_wrapper(SEXP triangle_in, SEXP n_boot, SEXP resids_type, SEXP boot_type, SEXP dist) {

    SEXP dim = getAttrib(triangle_in, R_DimSymbol);
    int n_dev = INTEGER(dim)[0];

    SEXP reserve;

    PROTECT(reserve = allocVector(REALSXP, asInteger(n_boot)));

    double triangle[n_dev][n_dev];

    for (int i=0; i<n_dev; i++) {
      for (int j=0; j<n_dev; j++) {
        triangle[j][i] = REAL(triangle_in)[j*n_dev + i];
      }
    }

    F77_CALL(mack_boot)(n_dev, triangle, asInteger(resids_type), asInteger(boot_type), asInteger(dist), asInteger(n_boot), REAL(reserve));

    UNPROTECT(1);

    return(reserve);
};



void F77_NAME(glm_boot)(int n_dev, double triangle[n_dev][n_dev], int n_boot, double reserve[n_boot]);

SEXP glm_boot_wrapper(SEXP triangle_in, SEXP n_boot) {

  SEXP dim = getAttrib(triangle_in, R_DimSymbol);
  int n_dev = INTEGER(dim)[0];

  SEXP reserve;

  PROTECT(reserve = allocVector(REALSXP, asInteger(n_boot)));

  double triangle[n_dev][n_dev];

  for (int i=0; i<n_dev; i++) {
    for (int j=0; j<n_dev; j++) {
      triangle[j][i] = REAL(triangle_in)[j*n_dev + i];
    }
  }

  F77_CALL(glm_boot)(n_dev, triangle, asInteger(n_boot), REAL(reserve));

  UNPROTECT(1);

  return(reserve);
};

F77_NAME();

SEXP mack_sim_wrapper(SEXP triangle, SEXP nboot, SEXP config, SEXP type) {

  SEXP dim = Rf_getAttrib(triangle, R_DimSymbol);
  int n_dev = INTEGER(dim)[0];



}

static const R_CallMethodDef callRoutines[] = {
    {"mack_boot_wrapper", (DL_FUNC) &mack_boot_wrapper, 5},
    {"glm_boot_wrapper", (DL_FUNC) &glm_boot_wrapper, 2},
    {NULL, NULL, 0}
    };

void R_init_patternBreak(DllInfo* dll) {
  R_registerRoutines(dll, NULL, callRoutines, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
};