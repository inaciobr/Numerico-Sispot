#include "node.h"

int leituraNode(node **nodes, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return 1;
    }

    int N;
    fscanf(fp, "%d", &N);

    node *nodesArquivo = malloc(N * sizeof(node));
    if (nodesArquivo == NULL) {
        printf("Não foi possível alocar memória no sistema.");
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
