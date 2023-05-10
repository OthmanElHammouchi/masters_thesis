#ifndef CONSTANTS
#define CONSTANTS

#include <Rcpp.h>
#include <map>

inline const int NONE = 0;

inline const int NORMAL = 1;
inline const int GAMMA = 2;

inline const int PARAMETRIC = 1;
inline const int RESID = 2;
inline const int PAIRS = 3;

inline const int STANDARDISED = 1;
inline const int MODIFIED = 2;
inline const int STUDENTISED = 3;
inline const int LOGNORMAL = 4;

inline const int SINGLE = 1;
inline const int CALENDAR = 2;
inline const int ORIGIN = 3;

inline Rcpp::List key = Rcpp::List::create(
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