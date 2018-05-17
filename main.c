#include <stdio.h>

#include "rede.h"
#include "matriz.h"
#include "numerico.h"

int main() {
    rede redePotencia;
    redePotencia = leituraRede("2_Reticulada");



    fluxoDePotenciaNewton(redePotencia);

    for (int k = 0; k < redePotencia.numBarras; k++) {
        printf("%d: %2.6f, %2.4f\n", redePotencia.barras[k].id, redePotencia.barras[k].tensao/redePotencia.barras[k].tensaoNominal, rad2Graus(redePotencia.barras[k].anguloTensao));
    }

    //testesZeroNewton();

    return 0;
}
