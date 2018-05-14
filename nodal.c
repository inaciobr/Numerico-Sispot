#include "nodal.h"

void leituraNodal(nodal *mNodal, int tamanhoMatriz, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return;
    }

    int N;
    fscanf(fp, "%d", &N);

    mNodal->condutancia = criaMatriz(tamanhoMatriz, tamanhoMatriz);
    mNodal->susceptancia = criaMatriz(tamanhoMatriz, tamanhoMatriz);

    int i, j;
    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d %d", &i, &j);
        fscanf(fp, "%lf %lf", &mNodal->condutancia.elemento[i][j],
                              &mNodal->susceptancia.elemento[i][j]);
    }

    fclose(fp);
}

void freeNodal(nodal *mNodal) {
    freeMatriz(&mNodal->condutancia);
    freeMatriz(&mNodal->susceptancia);
}
