#ifndef NUMERICO_H_INCLUDED
#define NUMERICO_H_INCLUDED

#include <math.h>

#include "matriz.h"

<<<<<<< HEAD
#define PI 3.14159265358979323846
#define MAX_ITERACOES 10
#define TOLERANCIA 1E-6

=======
>>>>>>> 6199045bbb7a8032d11bea93262d1e5dbeb1b311
matriz matrizFuncao(int numX, double x[], int numF, double* (*F[])(double[]));
void zeroNewton(int numX, double x[], double* (*F)(double[]), int numF, double* (*dF[])(double[]));

#endif // NUMERICO_H_INCLUDED
