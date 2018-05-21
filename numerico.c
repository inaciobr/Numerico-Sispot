#include "numerico.h"

/**
 * Monta a matriz coluna de termo conhecido do m�todo de newton.
 * Recebe como argumento a matriz M onde os resultados ser�o salvos
 * O vetor x de argumentos da fun��o utilizada
 * E um ponteiro para a fun��o F que ser� utilizada para montar a matriz.
 */
matriz matrizFuncao(matriz M, double x[], void (*F)(double[], double[])) {
    double *valores = calloc(M.numLinhas, sizeof(double));
    F(x, valores);

    for (int i = 0; i < M.numLinhas; i++)
        M.elemento[i][0] = valores[i];

    free(valores);

    return M;
}

/**
 * Monta a matriz Jacobiana do m�todo de newton.
 * Recebe como argumento a matriz M onde os resultados ser�o salvos
 * O vetor x de argumentos da fun��o utilizada
 * E um vetor de ponteiros para as fun��es F que ser�o utilizada para montar a matriz.
 */
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

/**
 * M�todo de newton para zero de fun��es.
 * Recebe como argumento o n�mero de vari�veis utilizadas,
 * um vetor x de valores iniciais e onde ser�o setados os valores
 * finais do resultado do m�todo de newton,
 * um ponteiro para uma fun��o F que se deseja zerar
 * e um vetor de ponteiros para as derivadas da fun��o F.
 */
int zeroNewton(int numX, double x[], void (*F)(double[], double[]), int numF, void (*dF[])(double[], double[])) {
    matriz Jx = criaMatriz(numF, numX);
    matriz Fx = criaMatriz(numX, 1);
    matriz R;

    int k;
    for (k = 0; k < MAX_ITERACOES; k++) {
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

    return k;
}

/**
 * Verifica a precis�o do m�todo de newton de forma relativa pela constante TOLERANCIA.
 * Se o valor absoluto de algum termo no vetor x multiplicado pela tolerancia definida
 * for maior do que a varia��o em uma itera��o do m�todo de newton, as itera��es continuam,
 * caso contr�rio o m�todo de newton � finalizado.
 */
int tolerancia(matriz R, double x[]) {
    for (int i = 0; i < R.numLinhas; i++)
        if (fabs(R.elemento[i][0]) > TOLERANCIA * fabs(x[i]))
            return 0;

    return 1;
}

/**
 * Converte um �ngulo de radianos para graus.
 */
double rad2Graus(double angulo) {
    return angulo * 180 / PI;
}
