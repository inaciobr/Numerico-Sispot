#include "numerico.h"

/**
 * M�todo de newton para zero de fun��es.
 * Recebe como argumento o n�mero de vari�veis utilizadas,
 * um vetor x de valores iniciais e onde ser�o setados os valores
 * finais do resultado do m�todo de newton,
 * um ponteiro para uma fun��o F que se deseja zerar
 * e um ponteiro para uma fun��o J referente � matriz jacobiana.
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

		/* Verifica a toler�ncia do resultado */
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
