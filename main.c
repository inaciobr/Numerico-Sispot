#include <stdio.h>

#include "rede.h"
#include "matriz.h"
#include "numerico.h"

int main() {
    rede redePotencia;
    redePotencia = leituraRede("1_Stevenson");

    double *x = calloc(2*redePotencia.numPQ + redePotencia.numPV, sizeof(double));

    int jMatriz = 0;
    for (int k = 0; k < redePotencia.numBarras; k++) {
        if (redePotencia.barras[k].tipo != B_PQ)
            continue;

        x[redePotencia.numPQ + redePotencia.numPV + jMatriz] = redePotencia.barras[k].tensao;
        x[jMatriz] = redePotencia.barras[k].anguloTensao;

        jMatriz++;
    }

    jMatriz = 0;
    for (int k = 0; k < redePotencia.numBarras; k++) {
        if (redePotencia.barras[k].tipo != B_PV)
            continue;

        x[redePotencia.numPQ + jMatriz] = redePotencia.barras[k].anguloTensao;

        jMatriz++;
    }

    for (int k = 0; k < 1; k++) {
        matriz Fx = funcaoDesvio(x, redePotencia);
        matriz Jx = jacobianaDesvios(redePotencia);

        matriz R = resolveSistemaLinear(Jx, Fx);

        for (int i = 0; i < R.numLinhas; i++)
            x[i] += R.elemento[i][0];
    }

    for (int k = 0; k < 2*redePotencia.numPQ + redePotencia.numPV; k++) {
        printf("%f\n", x[k]);
    }


    //testesZeroNewton();

    return 0;
}
