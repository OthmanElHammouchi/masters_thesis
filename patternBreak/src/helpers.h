#ifndef HELPERS
#define HELPERS

#include <stdlib.h>

inline double** array_create(int n, int m) {
    double* data = (double*)malloc(n * m * sizeof(double));
    double** array = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) {
      array[i] = &(data[m * i]);
    }
    return(array);
}

#endif