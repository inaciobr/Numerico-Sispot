#include "numerico.h"

matriz matrizFuncao(int numX, double x[], int numF, double *(*F[])(double[])) {
    matriz funcao = criaMatriz(numF, numX);

    for (int i = 0; i < numF; i++)
        for (int j = 0; j < numX; j++)
            funcao.elemento[i][j] = (F[i](x))[j];

    return funcao;
}

void zeroNewton(int numX, double x[], double *(*F)(double[]), int numF, double *(*dF[])(double[])) {

    for (int k = 0; k < 5; k++) {
        matriz Fx = transposta(matrizFuncao(numX, x, 1, &F));
        matriz Jx = matrizFuncao(numX, x, numF, dF);

        matriz R = produtoMatriz(inversa(Jx), multiplicaConstante(Fx, -1));


        for (int i = 0; i < R.numLinhas; i++) {
                x[i] += R.elemento[i][0];
        }
    }


    for (int i = 0; i < numX; i++) {
            printf("%15.8lf ", x[i]);

        printf("\n");
    }
}
