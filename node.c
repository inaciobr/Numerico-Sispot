#include "node.h"

int leituraNode(node **nodes, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("N�o foi poss�vel ler o arquivo %s.", arquivo);
        return 1;
    }

    int N;
    fscanf(fp, "%d", &N);

    node *nodesArquivo = malloc(N * sizeof(node));
    if (nodesArquivo == NULL) {
        printf("N�o foi poss�vel alocar mem�ria no sistema.");
        return 2;
    }

    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d %d %lf %lf",
                &nodesArquivo[i].posJ,
                &nodesArquivo[i].posK,
                &nodesArquivo[i].condutancia,
                &nodesArquivo[i].susceptancia);
    }

    fclose(fp);
    *nodes = nodesArquivo;
    return N;
}
