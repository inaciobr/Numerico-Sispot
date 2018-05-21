#include "numerico.h"

matriz matrizFuncao(matriz M, double x[], void (*F)(double[], double[])) {
    double *valores = calloc(M.numLinhas, sizeof(double));
    F(x, valores);

    for (int i = 0; i < M.numLinhas; i++)
        M.elemento[i][0] = valores[i];

    free(valores);

    return M;
}

matriz jacobiana(matriz M, double x[], void (*F[])(double[], double[])) {
    double *valores = calloc(M.numColunas, sizeof(double));

    for (int i = 0; i < M.numLinhas; i++) {
        F[i](x, valores);

        for (int j = 0; j < M.numColunas; j++)
            M.elemento[i][j] = valores[j];
    }

    free(valores);

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

        if (tolerancia(R, x)) {
            freeMatriz(&R);
            break;
        }

        freeMatriz(&R);
    }

    freeMatriz(&Jx);
    freeMatriz(&Fx);
}

int tolerancia(matriz R, double x[]) {
    double absX = 0.0, absR = 0.0;

    for (int i = 0; i < R.numLinhas; i++) {
        absR = R.elemento[i][0] > 0 ? R.elemento[i][0] : -R.elemento[i][0];
        absX = x[i] > 0 ? x[i] : -x[i];

        if (TOLERANCIA * absX > absR)
            return 1;
    }

    return 0;
}

double rad2Graus(double angulo) {
    return angulo * 180 / PI;
}

// TESTES INICIAIS
void F1x(double x[2], double r[2]) {
    r[0] = 2*(x[0] - 2);
    r[1] = 2*(x[1] - 3);
}

void dF1dx1(double x[2], double r[2]) {
    r[0] = 2;
    r[1] = 0;
}

void dF1dx2(double x[2], double r[2]) {
    r[0] = 0;
    r[1] = 2;
}

void testesZeroNewton() {
    double x1[2] = {0.1, -0.1};
    void (*dF1[2]) = {&dF1dx1, &dF1dx2};

    zeroNewton(2, x1, &F1x, 2, &dF1);

    for (int i = 0; i < 2; i++)
        printf("%15.10f ", x1[i]);

    printf("\n\n");

}
