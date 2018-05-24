#include "numerico.h"

/**
 * Método de newton para zero de funções.
 * Recebe como argumento o número de variáveis utilizadas,
 * um vetor x de valores iniciais e onde serão setados os valores
 * finais do resultado do método de newton,
 * um ponteiro para uma função F que se deseja zerar
 * e um ponteiro para uma função J referente à matriz jacobiana.
 */
int zeroNewton(int numX, double x[], void (*F)(matriz*, double[]), void (*J)(matriz*, double[])) {
    matriz Jx = criaMatriz(numX, numX);
    matriz Fx = criaMatriz(numX, 1);
    matriz R;

    int k;
    for (k = 0; k < 5; k++) {
        F(&Fx, x);
        J(&Jx, x);

        R = resolveSistemaLinear(Jx, multiplicaConstante(Fx, -1));

		/* Atualiza o vetor x */
        for (int i = 0; i < R.numLinhas; i++)
            x[i] += R.elemento[i][0];

		/* Verifica a tolerância do resultado */
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
