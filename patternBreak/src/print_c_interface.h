#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rinterface.h>

// C wrappers for R's print functionality.

void F77_SUB(Rprintf)(char* string) {
    Rprintf(string);
};

void F77_SUB(R_Flushconsole)(char* string) {
  R_FlushConsole();
};