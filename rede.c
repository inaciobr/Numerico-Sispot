#include "rede.h"

rede leituraRede(char arquivo[]) {
    rede r;

    r.numPQ = 0;
    r.numPV = 0;
    r.numSwing = 0;

    char caminho1[256], caminho2[256];

    sprintf(caminho1, "Redes/%s/%s_DadosBarras.txt", arquivo, arquivo);
    sprintf(caminho2, "Redes/%s/%s_YNodal.txt", arquivo, arquivo);

    leituraBarra(&r, caminho1);
    leituraNodal(&r, caminho2);

    return r;
}

void leituraBarra(rede *r, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return;
    }

    fscanf(fp, "%d", &r->numBarras);

    barra *barrasArquivo = malloc(r->numBarras * sizeof(barra));

    for (int i = 0; i < r->numBarras; i++) {
        fscanf(fp, "%d %d %lf", &barrasArquivo[i].id,
                                &barrasArquivo[i].tipo,
                                &barrasArquivo[i].tensaoNominal);

        /* Valores especificados para a rede em cada caso */
        switch (barrasArquivo[i].tipo) {
        case B_PQ:
            fscanf(fp, "%lf %lf", &barrasArquivo[i].potenciaAtiva,
                                  &barrasArquivo[i].potenciaReativa);

            barrasArquivo[i].anguloTensao = 0.0;
            barrasArquivo[i].tensao = barrasArquivo[i].tensaoNominal;
            barrasArquivo[i].potenciaAtiva = -barrasArquivo[i].potenciaAtiva;
            r->numPQ++;

            break;

        case B_PV:
            fscanf(fp, "%lf %lf", &barrasArquivo[i].potenciaAtiva,
                                  &barrasArquivo[i].tensao);

            barrasArquivo[i].anguloTensao = 0.0;
            r->numPV++;

            break;

        case B_SWING:
            fscanf(fp, "%lf %lf", &barrasArquivo[i].tensao,
                                  &barrasArquivo[i].anguloTensao);

            barrasArquivo[i].potenciaAtiva = 0.0;
            barrasArquivo[i].potenciaReativa = 0.0;
            r->numSwing++;

            break;
        }
    }

    fclose(fp);
    r->barras = barrasArquivo;
}

void leituraNodal(rede *r, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("Não foi possível ler o arquivo %s.", arquivo);
        return;
    }

    int N;
    fscanf(fp, "%d", &N);

    nodal mNodal;
    mNodal.condutancia = criaMatriz(r->numBarras, r->numBarras);
    mNodal.susceptancia = criaMatriz(r->numBarras, r->numBarras);

    int linha, coluna;
    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d %d", &linha, &coluna);
        fscanf(fp, "%lf %lf", &mNodal.condutancia.elemento[linha][coluna],
                              &mNodal.susceptancia.elemento[linha][coluna]);
    }

    r->mNodal = mNodal;

    fclose(fp);
}

void freeNodal(nodal *mNodal) {
    freeMatriz(&mNodal->condutancia);
    freeMatriz(&mNodal->susceptancia);
}


void fP(double resultado[], rede r) {
    int iFP = 0;
    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo == B_SWING)
            continue;

        resultado[iFP] = 0.0;
        for (int k = 0; k < r.numBarras; k++) {
            double angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;

            resultado[iFP] += r.barras[k].tensao * (r.mNodal.condutancia.elemento[j][k]*cos(angulo)
                                                  - r.mNodal.susceptancia.elemento[j][k]*sin(angulo));
        }

        resultado[iFP] *= r.barras[j].tensao;

        if (r.barras[j].tipo == B_PV)
            resultado[iFP] -= r.barras[j].potenciaAtiva;

        iFP++;
    }
}

void fQ(double resultado[], rede r) {
    int iFP = 0;
    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo != B_PQ)
            continue;

        resultado[iFP] = 0.0;
        for (int k = 0; k < r.numBarras; k++) {
            double angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;

            resultado[iFP] += r.barras[k].tensao * (r.mNodal.condutancia.elemento[j][k]*cos(angulo)
                                                  - r.mNodal.susceptancia.elemento[j][k]*sin(angulo));
        }

        resultado[iFP] *= r.barras[j].tensao;

        iFP++;
    }
}

matriz funcaoDesvio(double x[], rede r) {
    matriz M = criaMatriz(2*r.numPQ + r.numPV, 1);

    double *theta = x;
    double *tensao = &x[r.numPQ + r.numPV];

    int j = 0, k = r.numPQ;
    for (int i = 0; i < r.numBarras; i++) {
        if (r.barras[i].tipo == B_PQ) {
            r.barras[i].tensao = tensao[j];
            r.barras[i].anguloTensao = theta[j];
            j++;
        } else if (r.barras[i].tipo == B_PV) {
            r.barras[i].anguloTensao = theta[k];
            k++;
        }
    }

    double *resultado = calloc(2*r.numPQ + r.numPV, sizeof(double));

    fP(resultado, r);
    fQ(&resultado[r.numPQ + r.numPV], r);

    for (int i = 0; i < 2*r.numPQ + r.numPV; i++)
        M.elemento[i][0] = -resultado[i];

    return M;
}


//double *theta = x;
//double *tensao = &x[r.numPQ + r.numPV];
matriz jacobianaDesvios(rede r) {
    matriz M = criaMatriz(2*r.numPQ + r.numPV, 2*r.numPQ + r.numPV);
    double angulo;

    //del FP / del THETA => (r.numPQ + r.numPV) x (r.numPQ + r.numPV)
    int jMatriz = 0, iFP = 0;
    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo == B_SWING)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r.numBarras; k++) {
            if (r.barras[k].tipo == B_SWING)
                continue;

            if (iFP == jMatriz) {
                M.elemento[iFP][jMatriz] = 0.0;

                for (int i = 0; i < r.numBarras; i++) {
                    if(i == j)
                        continue;
                    angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                    M.elemento[iFP][jMatriz] += r.barras[i].tensao
                                                * (r.mNodal.condutancia.elemento[j][i]*sin(angulo)
                                                    + r.mNodal.susceptancia.elemento[j][i]*cos(angulo));
                }

                M.elemento[iFP][jMatriz] *= r.barras[j].tensao;
            } else {
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[iFP][jMatriz] = -r.barras[j].tensao * r.barras[k].tensao
                                            * (r.mNodal.condutancia.elemento[j][k]*sin(angulo)
                                                + r.mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }

    //del FP / del V => (r.numPQ + r.numPV) x (r.numPQ)
    iFP = 0;
    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo == B_SWING)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r.numBarras; k++) {
            if (r.barras[k].tipo != B_PQ)
                continue;

            if (iFP == jMatriz) {

                M.elemento[iFP][r.numPQ + r.numPV + jMatriz] = 0.0;
                for (int i = 0; i < r.numBarras; i++) {
                    if(i == j){
                        M.elemento[iFP][r.numPQ + r.numPV + jMatriz] += 2*r.barras[j].tensao*r.mNodal.condutancia.elemento[j][j];
                    } else {
                        angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                        M.elemento[iFP][r.numPQ + r.numPV + jMatriz] += r.barras[i].tensao * (r.mNodal.condutancia.elemento[j][i]*cos(angulo) - r.mNodal.susceptancia.elemento[j][i]*sin(angulo));
                    }
                }

            } else {
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[iFP][r.numPQ + r.numPV + jMatriz] = r.barras[j].tensao * (r.mNodal.condutancia.elemento[j][k]*cos(angulo) - r.mNodal.susceptancia.elemento[j][k]*sin(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }


    //del FQ / del THETA => (r.numPQ) x (r.numPQ + r.numPV)
    iFP = 0;
    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo != B_PQ)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r.numBarras; k++) {
            if (r.barras[k].tipo == B_SWING)
                continue;

            if (iFP == jMatriz) {
                M.elemento[r.numPQ + r.numPV + iFP][jMatriz] = 0.0;

                for (int i = 0; i < r.numBarras; i++) {
                    if(i == j)
                        continue;

                    angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                    M.elemento[r.numPQ + r.numPV + iFP][jMatriz] += r.barras[i].tensao*(r.mNodal.condutancia.elemento[j][i]*cos(angulo) - r.mNodal.susceptancia.elemento[j][i]*sin(angulo));
                }

                M.elemento[r.numPQ + r.numPV + iFP][jMatriz] *= r.barras[j].tensao;
            } else {
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[r.numPQ + r.numPV + iFP][jMatriz] = -r.barras[j].tensao * r.barras[k].tensao
                                                            * (r.mNodal.condutancia.elemento[j][k]*cos(angulo)
                                                                - r.mNodal.susceptancia.elemento[j][k]*sin(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }


    //del FQ / del V => (r.numPQ) x (r.numPQ)
    iFP = 0;
    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo != B_PQ)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r.numBarras; k++) {
            if (r.barras[k].tipo != B_PQ)
                continue;

            if (iFP == jMatriz) {
                M.elemento[r.numPQ + r.numPV + iFP][jMatriz] = 0.0;

                for (int i = 0; i < r.numBarras; i++) {
                    if(i == j){
                        M.elemento[r.numPQ + r.numPV + iFP][r.numPQ + r.numPV + jMatriz] += -2*r.barras[j].tensao*r.mNodal.susceptancia.elemento[j][j];
                    } else {
                    angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                    M.elemento[r.numPQ + r.numPV + iFP][r.numPQ + r.numPV + jMatriz] += -r.barras[i].tensao
                                                                                    * (r.mNodal.condutancia.elemento[j][i]*sin(angulo)
                                                                                        + r.mNodal.susceptancia.elemento[j][i]*cos(angulo));
                     }
                }
            } else {
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[r.numPQ + r.numPV + iFP][r.numPQ + r.numPV + jMatriz] = -r.barras[j].tensao
                                                                                * (r.mNodal.condutancia.elemento[j][k]*sin(angulo)
                                                                                    + r.mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }

    return M;

}
