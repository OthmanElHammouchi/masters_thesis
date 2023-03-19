#ifndef CONSTANTS
#define CONSTANTS

#include <Rcpp.h>
#include <unordered_map>

const int NONE = 0;

const int NORMAL = 1;
const int GAMMA = 2;

const int PARAMETRIC = 1;
const int NON_PARAMETRIC = 2;
const int PAIRS = 3;

const int STANDARDISED = 1;
const int MODIFIED = 2;
const int STUDENTISED = 3;
const int LOGNORMAL = 4;

const int SINGLE = 1;
const int CALENDAR = 2;
const int ORIGIN = 3;

std::unordered_map<Rcpp::String, int> key = {
  {"none", NONE},
  {"normal", NORMAL},
  {"gamma", GAMMA},
  {"parametric", PARAMETRIC},
  {"non-parametric", NON_PARAMETRIC},
  {"pairs", PAIRS},
  {"standardised", STANDARDISED},
  {"modified", MODIFIED},
  {"studentised", STUDENTISED},
  {"log-normal", LOGNORMAL},
  {"single", SINGLE},
  {"calendar", CALENDAR},
  {"origin", ORIGIN}
};

#endif