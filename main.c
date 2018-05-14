#include <stdio.h>

#include "node.h"
#include "barra.h"
#include "matriz.h"
#include "numerico.h"

int main() {
    node nodes;
    barra barras;

    int numBarras = leituraBarra(&barras, "Redes/1_Stevenson/1_Stevenson_DadosBarras.txt");
    leituraNode(&nodes, numBarras, "Redes/1_Stevenson/1_Stevenson_Ynodal.txt");

    for (int i = 0; i < numBarras; i++) {
        for (int j = 0; j < numBarras; j++) {
            printf("%f ", nodes.condutancia[i][j]);
        }
        printf("\n");
    }


    return 0;
}
