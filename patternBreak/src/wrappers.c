#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rinterface.h>

#include "rng_c_interface.h"

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

void F77_NAME(mack_sim)(int n_dev, double triangle[n_dev][n_dev], int n_config, int m_config, double config[n_config][m_config], int type, int n_boot, double results[n_config * n_boot][m_config + 1]);

SEXP mack_sim_wrapper(SEXP triangle_in, SEXP n_boot, SEXP config_in, SEXP type) {

  SEXP dim = Rf_getAttrib(triangle_in, R_DimSymbol);
  int n_dev = INTEGER(dim)[0];

  dim = Rf_getAttrib(config_in, R_DimSymbol);
  int n_config = INTEGER(dim)[0];
  int m_config = INTEGER(dim)[1];

  double triangle[n_dev][n_dev];

  for (int i=0; i<n_dev; i++) {
    for (int j=0; j<n_dev; j++) {
      triangle[j][i] = REAL(triangle_in)[j*n_dev + i];
    }
  }

  double (*config)[n_config] = malloc(sizeof(double[m_config][n_config]));

  for (int i=0; i<n_config; i++) {
    for (int j=0; j<m_config; j++) {
      config[j][i] = REAL(config_in)[j*n_config + i];
    }
  }

  int n_results = asInteger(n_boot) * n_config;
  int m_results = m_config + 1;

  double (*results_in)[n_results] = malloc(sizeof(double[m_results][n_results]));

  F77_CALL(mack_sim)(n_dev, triangle, n_config, m_config, config, asInteger(type), asInteger(n_boot), results_in);

  SEXP results;

  PROTECT(results = allocVector(REALSXP, m_results * n_results));

    for (int i=0; i<n_results; i++) {
      for (int j=0; j<m_results; j++) {
        REAL(results)[j*n_results + i] = results_in[j][i];
      }
    }

  SEXP results_dim = PROTECT(allocVector(INTSXP, 2));
  INTEGER(results_dim)[0] = n_results;
  INTEGER(results_dim)[1] = m_results;

  Rf_setAttrib(results, R_DimSymbol, results_dim);

  UNPROTECT(2);

  free(config);
  free(results_in);

  return(results);
};

static const R_CallMethodDef callRoutines[] = {
    {"mack_boot_wrapper", (DL_FUNC) &mack_boot_wrapper, 5},
    {"glm_boot_wrapper", (DL_FUNC) &glm_boot_wrapper, 2},
    {"mack_sim_wrapper", (DL_FUNC) &mack_sim_wrapper, 4},
    {NULL, NULL, 0}
    };

void R_init_patternBreak(DllInfo* dll) {
  R_registerRoutines(dll, NULL, callRoutines, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
};