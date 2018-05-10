#include "numerico.h"

matriz matrizJacobiana(double *x, double* (*dF1)(double x[]), double* (*dF2)(double x[]), double* (*dF3)(double x[])) {
    matriz Jacobiana = criaMatriz(3, 3);

    Jacobiana.elemento[0] = dF1(x);
    Jacobiana.elemento[1] = dF2(x);
    Jacobiana.elemento[2] = dF3(x);

    return Jacobiana;
}
