#include "barra.h"

int leituraBarra(barra **barras, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("N�o foi poss�vel ler o arquivo %s.", arquivo);
        return 1;
    }

    int N;
    fscanf(fp, "%d", &N);

    barra *barrasArquivo = malloc(N * sizeof(barra));
    if (barrasArquivo == NULL) {
        printf("N�o foi poss�vel alocar mem�ria no sistema.");
        return 2;
    }

    for (int i = 0; i < N; i++) {

    }

    fclose(fp);
    *barras = barrasArquivo;
    return N;
}
