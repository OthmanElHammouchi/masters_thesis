#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rinterface.h>

// C wrappers for R's RNG functionality

void F77_SUB(GetRNGstate_fentry)(void) { GetRNGstate(); };

void F77_SUB(PutRNGstate_fentry)(void) { PutRNGstate(); };

double F77_SUB(rgamma_fentry)(double shape, double scale) {
  return(rgamma(shape, scale));
};

double F77_SUB(rnorm_fentry)(double mean, double sd) {
  return(rnorm(mean, sd));
};

double F77_SUB(rpois_fentry)(double lambda) {
  return(rpois(lambda));
};