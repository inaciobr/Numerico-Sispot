#include "numerico.h"

matriz matrizFuncao(int numX, double x[], int numF, ...) {
    va_list valist;
    matriz funcao = criaMatriz(numF, numX);

    va_start(valist, numF);

    for (int i = 0; i < numF; i++) {
        double *(*F)(double[]) = va_arg(valist, double* (*));
        double *v = F(x);

        for (int j = 0; j < numX; j++)
            funcao.elemento[i][j] = v[j];

        free(v);
    }

    va_end(valist);

    return funcao;
}

double* zeroNewton(int numX, double x[], int numF, ...) {

}
