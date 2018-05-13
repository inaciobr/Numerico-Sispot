#include "numerico.h"

matriz matrizFuncao(matriz M, double x[], void (*F)(double[], double[])) {
    double *valores = malloc(M.numLinhas * sizeof(double));

    for (int i = 0; i < M.numColunas; i++) {
        F(x, valores);
        M.elemento[i][0] = valores[i];
    }

    return M;
}

matriz jacobiana(matriz M, double x[], void (*F[])(double[], double[])) {
    double *valores = malloc(M.numColunas * sizeof(double));

    for (int i = 0; i < M.numLinhas; i++) {
        F[i](x, valores);

        for (int j = 0; j < M.numColunas; j++)
            M.elemento[i][j] = valores[j];
    }

    return M;
}

void zeroNewton(int numX, double x[], void (*F)(double[], double[]), int numF, void (*dF[])(double[], double[])) {
    matriz Jx = criaMatriz(numF, numX);
    matriz Fx = criaMatriz(numX, 1);
    matriz R;

    for (int k = 0; k < MAX_ITERACOES; k++) {
        Fx = matrizFuncao(Fx, x, F);
        Jx = jacobiana(Jx, x, dF);
        R = resolveSistemaLinear(Jx, multiplicaConstante(Fx, -1));

        for (int i = 0; i < R.numLinhas; i++)
                x[i] += R.elemento[i][0];
    }

    freeMatriz(&Jx);
    freeMatriz(&Fx);
}
