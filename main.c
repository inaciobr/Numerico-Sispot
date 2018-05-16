#include <stdio.h>

#include "rede.h"
#include "matriz.h"
#include "numerico.h"

int main() {
    rede redePotencia;
    redePotencia = leituraRede("1_Stevenson");
    matriz J = jacobianaDesvios(redePotencia);
    //double res[2*redePotencia.numPQ + redePotencia.numPV];

    // De 0 a redePotencia.numPQ + redePotencia.numPV => theta, após => tensão
    //double *x = calloc(2*redePotencia.numPQ + redePotencia.numPV, sizeof(double));

    //funcaoDesvio(x, res, redePotencia);
    //matriz M = jacobianaDesvios(x, res, redePotencia);

    //for (int i = 0; i < 2*redePotencia.numPQ + redePotencia.numPV; i++)
      //  printf("%f ", res[i]);
    printMatriz(J);
    return 0;
}
