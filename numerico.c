#include "numerico.h"

matriz matrizFuncao(int numX, double x[], int numF, double* (*F[])(double[])) {
    matriz funcao = criaMatriz(numF, numX);
    double *valores;

    for (int i = 0; i < numF; i++) {
        valores = F[i](x);

        for (int j = 0; j < numX; j++)
            funcao.elemento[i][j] = valores[j];
    }

    return funcao;
}

void zeroNewton(int numX, double x[], double* (*F)(double[]), int numF, double *(*dF[])(double[])) {
    matriz Fx, Jx, R;

    for (int k = 0; k < MAX_ITERACOES; k++) {
        Fx = transposta(matrizFuncao(numX, x, 1, &F));
        Jx = matrizFuncao(numX, x, numF, dF);

        R = resolveSistemaLinear(Jx, multiplicaConstante(Fx, -1));


        for (int i = 0; i < R.numLinhas; i++) {
                x[i] += R.elemento[i][0];
        }

        freeMatriz(&Fx);
        freeMatriz(&Jx);
    }
}
