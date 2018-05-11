#ifndef NUMERICO_H_INCLUDED
#define NUMERICO_H_INCLUDED

#include <math.h>
#include <stdarg.h>

#include "matriz.h"

matriz matrizFuncao(int numX, double x[], int numF, ...);
double* zeroNewton(int numX, double x[], int numF, ...);

#endif // NUMERICO_H_INCLUDED
