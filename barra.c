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
        fscanf(fp, "%d %d %lf",
                &barrasArquivo[i].id,
                &barrasArquivo[i].tipo,
                &barrasArquivo[i].tensaoNominal);

        switch (barrasArquivo[i].tipo) {
        case B_PQ:
            fscanf(fp, "%lf %lf",
                    &barrasArquivo[i].potenciaAtiva,
                    &barrasArquivo[i].potenciaReativa);
            break;

        case B_PV:
            fscanf(fp, "%lf %lf",
                    &barrasArquivo[i].potenciaAtiva,
                    &barrasArquivo[i].tensaoNominal);
            barrasArquivo[i].potenciaAtiva = -barrasArquivo[i].potenciaAtiva;
            break;

        case B_SWING:
            fscanf(fp, "%lf %lf",
                    &barrasArquivo[i].tensaoNominal,
                    &barrasArquivo[i].anguloTensao);
            break;
        }
    }

    fclose(fp);
    *barras = barrasArquivo;
    return N;
}
