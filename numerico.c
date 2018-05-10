#include "numerico.h"

matriz matrizFuncao(int numX, double *x, int numF, ...) {
    va_list valist;
    matriz funcao = criaMatriz(numF, numX);

    va_start(valist, numF);

    for (int i = 0; i < numF; i++) {
        double* (*F)(double[]) = va_arg(valist, double* (*));
        funcao.elemento[i] = F(x);
    }

    va_end(valist);

    return funcao;
}
