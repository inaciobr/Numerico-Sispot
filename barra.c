#include "barra.h"

int leituraBarra(char arquivo[]) {
    FILE *fp;

    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s", arquivo);
        return 1;
    }

    fclose(fp);
    return 0;
}
