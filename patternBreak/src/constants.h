#ifndef CONSTANTS
#define CONSTANTS

#include <Rcpp.h>
#include <map>

const int NONE = 0;

const int NORMAL = 1;
const int GAMMA = 2;

const int PARAMETRIC = 1;
const int RESID = 2;
const int PAIRS = 3;

const int STANDARDISED = 1;
const int MODIFIED = 2;
const int STUDENTISED = 3;
const int LOGNORMAL = 4;

const int SINGLE = 1;
const int CALENDAR = 2;
const int ORIGIN = 3;

Rcpp::List key = Rcpp::List::create(
  Rcpp::Named("") = NONE,
  Rcpp::Named("normal") = NORMAL,
  Rcpp::Named("gamma") = GAMMA,
  Rcpp::Named("parametric") = PARAMETRIC,
  Rcpp::Named("residuals") = RESID,
  Rcpp::Named("pairs") = PAIRS,
  Rcpp::Named("standardised") = STANDARDISED,
  Rcpp::Named("modified") = MODIFIED,
  Rcpp::Named("studentised") = STUDENTISED,
  Rcpp::Named("log-normal") = LOGNORMAL,
  Rcpp::Named("single") = SINGLE,
  Rcpp::Named("calendar") = CALENDAR,
  Rcpp::Named("origin") = ORIGIN
);

#endif