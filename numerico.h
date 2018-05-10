#ifndef NUMERICO_H_INCLUDED
#define NUMERICO_H_INCLUDED

#include <math.h>
#include <stdarg.h>

#include "matriz.h"

matriz matrizFuncao(int numX, double *x, int numF, ...);
matriz newton();

#endif // NUMERICO_H_INCLUDED
