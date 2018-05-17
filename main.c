#include <stdio.h>

#include "rede.h"
#include "matriz.h"
#include "numerico.h"

int main() {
    rede redePotencia;
    redePotencia = leituraRede("1_Stevenson");

    double *x = calloc(2*redePotencia.numPQ + redePotencia.numPV, sizeof(double));

    int jMatriz = 0, jPQ = 0;
    for (int k = 0; k < redePotencia.numBarras; k++) {
        if (redePotencia.barras[k].tipo == B_SWING)
            continue;

        x[jMatriz] = redePotencia.barras[k].anguloTensao;

        if (redePotencia.barras[k].tipo == B_PQ) {
            x[redePotencia.numPQ + redePotencia.numPV + jPQ] = redePotencia.barras[k].tensao;
            jPQ++;
        }

        jMatriz++;
    }

    for (int k = 0; k < 5; k++) {
        matriz Fx = funcaoDesvio(x, redePotencia);
        matriz Jx = jacobianaDesvios(redePotencia);


        matriz R = resolveSistemaLinear(Jx, Fx);

        for (int i = 0; i < R.numLinhas; i++)
            x[i] += R.elemento[i][0];
    }

    for (int k = 0; k < redePotencia.numBarras; k++) {
        printf("%d: %2.5f\n", redePotencia.barras[k].id, redePotencia.barras[k].tensao/redePotencia.barras[k].tensaoNominal);
    }

    //testesZeroNewton();

    return 0;
}
