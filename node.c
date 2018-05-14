#include "node.h"

void leituraNode(node *nodes, int tamanhoMatriz, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return;
    }

    int N;
    fscanf(fp, "%d", &N);

    node nodesArquivo;
    nodesArquivo.condutancia = malloc(tamanhoMatriz * sizeof(double *));
    nodesArquivo.susceptancia = malloc(tamanhoMatriz * sizeof(double *));

    for (int i = 0; i < tamanhoMatriz; i++) {
        nodesArquivo.condutancia[i] = calloc(tamanhoMatriz, sizeof(double));
        nodesArquivo.susceptancia[i] = calloc(tamanhoMatriz, sizeof(double));
    }

    int i, j;
    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d %d", &i, &j);
        fscanf(fp, "%lf %lf", &nodesArquivo.condutancia[i][j],
                              &nodesArquivo.susceptancia[i][j]);
    }

    fclose(fp);
    *nodes = nodesArquivo;
}

// freeNode;
