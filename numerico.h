#ifndef NUMERICO_H_INCLUDED
#define NUMERICO_H_INCLUDED

#include <math.h>

#include "matriz.h"

matriz matrizFuncao(int numX, double x[], int numF, double* (*F[])(double[]));
void zeroNewton(int numX, double x[], double* (*F)(double[]), int numF, double* (*dF[])(double[]));

#endif // NUMERICO_H_INCLUDED
