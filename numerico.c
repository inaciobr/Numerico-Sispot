#include "numerico.h"

matriz matrizFuncao(matriz M, double x[], double* (*F)(double[])) {
    double *valores;

    for (int i = 0; i < M.numColunas; i++) {
        valores = F(x);

        M.elemento[i][0] = valores[i];
    }

    return M;
}

matriz jacobiana(matriz M, double x[], double* (*F[])(double[])) {
    double *valores;

    for (int i = 0; i < M.numLinhas; i++) {
        valores = F[i](x);

        for (int j = 0; j < M.numColunas; j++)
            M.elemento[i][j] = valores[j];
    }

    return M;
}

void zeroNewton(int numX, double x[], double* (*F)(double[]), int numF, double *(*dF[])(double[])) {
    matriz Jx = criaMatriz(numF, numX);
    matriz Fx = criaMatriz(numX, 1);

    for (int k = 0; k < MAX_ITERACOES; k++) {
        Fx = matrizFuncao(Fx, x, F);
        Jx = jacobiana(Jx, x, dF);

        matriz R = resolveSistemaLinear(Jx, multiplicaConstante(Fx, -1));

        for (int i = 0; i < R.numLinhas; i++)
                x[i] += R.elemento[i][0];
    }

    freeMatriz(&Jx);
    freeMatriz(&Fx);
}
