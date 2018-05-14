#include "barra.h"

int leituraBarra(barra **barras, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return 0;
    }

    int N;
    fscanf(fp, "%d", &N);

    barra *barrasArquivo = malloc(N * sizeof(barra));

    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d %d %lf", &barrasArquivo[i].id,
                                &barrasArquivo[i].tipo,
                                &barrasArquivo[i].tensao);

        switch (barrasArquivo[i].tipo) {
        case B_PQ:
            fscanf(fp, "%lf %lf", &barrasArquivo[i].potenciaAtiva,
                                  &barrasArquivo[i].potenciaReativa);
            break;

        case B_PV:
            fscanf(fp, "%lf %lf", &barrasArquivo[i].potenciaAtiva,
                                  &barrasArquivo[i].tensao);
            barrasArquivo[i].potenciaAtiva = -barrasArquivo[i].potenciaAtiva;
            break;

        case B_SWING:
            fscanf(fp, "%lf %lf", &barrasArquivo[i].tensao,
                                  &barrasArquivo[i].anguloTensao);
            break;
        }
    }

    fclose(fp);
    *barras = barrasArquivo;
    return N;
}
