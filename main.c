#include <stdio.h>

#include "nodal.h"
#include "barra.h"
#include "matriz.h"
#include "numerico.h"

int main() {
    barra *barras;
    nodal mNodal;

    int numBarras = leituraBarra(&barras, "Redes/1_Stevenson/1_Stevenson_DadosBarras.txt");
    leituraNodal(&mNodal, numBarras, "Redes/1_Stevenson/1_Stevenson_Ynodal.txt");

    for (int i = 0; i < numBarras; i++) {
        for (int j = 0; j < numBarras; j++) {
            printf("%15.8f ", mNodal.condutancia.elemento[i][j]);
        }
        printf("\n");
    }

    return 0;
}
