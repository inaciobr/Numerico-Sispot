#ifndef NUMERICO_H_INCLUDED
#define NUMERICO_H_INCLUDED

#include <math.h>

#include "matriz.h"

matriz matrizJacobiana(double *x, double* (*dF1)(double x[]), double* (*dF2)(double x[]), double* (*dF3)(double x[]));

matriz matrizFuncao(matriz x, double* (*F)(double x[]));
matriz newton();

#endif // NUMERICO_H_INCLUDED
