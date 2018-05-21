#include "numerico.h"

/**
 * Monta a matriz coluna de termo conhecido do método de newton.
 * Recebe como argumento a matriz M onde os resultados serão salvos
 * O vetor x de argumentos da função utilizada
 * E um ponteiro para a função F que será utilizada para montar a matriz.
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
 * Monta a matriz Jacobiana do método de newton.
 * Recebe como argumento a matriz M onde os resultados serão salvos
 * O vetor x de argumentos da função utilizada
 * E um vetor de ponteiros para as funções F que serão utilizada para montar a matriz.
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
 * Método de newton para zero de funções.
 * Recebe como argumento o número de variáveis utilizadas,
 * um vetor x de valores iniciais e onde serão setados os valores
 * finais do resultado do método de newton,
 * um ponteiro para uma função F que se deseja zerar
 * e um vetor de ponteiros para as derivadas da função F.
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
 * Verifica a precisão do método de newton de forma relativa pela constante TOLERANCIA.
 * Se o valor absoluto de algum termo no vetor x multiplicado pela tolerancia definida
 * for maior do que a variação em uma iteração do método de newton, as iterações continuam,
 * caso contrário o método de newton é finalizado.
 */
int tolerancia(matriz R, double x[]) {
    for (int i = 0; i < R.numLinhas; i++)
        if (fabs(R.elemento[i][0]) > TOLERANCIA * fabs(x[i]))
            return 0;

    return 1;
}

/**
 * Converte um ângulo de radianos para graus.
 */
double rad2Graus(double angulo) {
    return angulo * 180 / PI;
}
