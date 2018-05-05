#include <stdio.h>
#include "node.h"
#include "barra.h"
#include "matriz.h"

int main() {
    node *nodes = NULL;
    barra *barras = NULL;
    int numNodes, numBarras;

    numNodes = leituraNode(&nodes, "Redes/1_Stevenson/1_Stevenson_Ynodal.txt");
    numBarras = leituraBarra(&barras, "Redes/1_Stevenson/1_Stevenson_DadosBarras.txt");

    printf("%d, %d", nodes[numNodes - 1].posJ, barras[numBarras - 1].id);

    free(nodes);
    nodes = NULL;
    free(barras);
    barras = NULL;

    return 0;
}
