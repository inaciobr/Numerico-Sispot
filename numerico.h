#ifndef NUMERICO_H_INCLUDED
#define NUMERICO_H_INCLUDED

#include <math.h>

#include "matriz.h"

#define PI 3.14159265358979323846
#define MAX_ITERACOES 10
#define TOLERANCIA 1E-6


matriz matrizFuncao(matriz M, double x[], void (*F)(double[], double[]));
matriz jacobiana(matriz M, double x[], void (*F[])(double[], double[]));
int zeroNewton(int numX, double x[], void (*F)(double[], double[]), int numF, void (*dF[])(double[], double[]));

int tolerancia(matriz R, double x[]);
double rad2Graus(double angulo);

#endif // NUMERICO_H_INCLUDED
