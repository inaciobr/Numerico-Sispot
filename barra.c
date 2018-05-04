#include "barra.h"

int leituraBarra(barra **barras, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return 1;
    }

    int N;
    fscanf(fp, "%d", &N);

    barra *barrasArquivo = malloc(N * sizeof(barra));
    if (barrasArquivo == NULL) {
        printf("Não foi possível alocar memória no sistema.");
        return 2;
    }

    for (int i = 0; i < N; i++) {

    }

    fclose(fp);
    *barras = barrasArquivo;
    return N;
}
