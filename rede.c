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

void fluxoDePotenciaNewton(rede r) {
    matriz Jx = criaMatriz(2*r.numPQ + r.numPV, 2*r.numPQ + r.numPV);
    matriz Fx = criaMatriz(2*r.numPQ + r.numPV, 1);
    matriz R;

    double *x = calloc(2*r.numPQ + r.numPV, sizeof(double));

    // Valores iniciais de x.
    int jMatriz = 0, jPQ = 0;
    for (int k = 0; k < r.numBarras; k++) {
        if (r.barras[k].tipo == B_SWING)
            continue;

        x[jMatriz] = r.barras[k].anguloTensao;

        if (r.barras[k].tipo == B_PQ) {
            x[r.numPQ + r.numPV + jPQ] = r.barras[k].tensao;
            jPQ++;
        }

        jMatriz++;
    }

    for (int k = 0; k < MAX_ITERACOES; k++) {
        Fx = funcaoDesvio(Fx, r);
        Jx = jacobianaDesvios(Jx, r);

        R = resolveSistemaLinear(Jx, Fx);

        for (int i = 0; i < R.numLinhas; i++)
            x[i] += R.elemento[i][0];

        atualizaRede(x, r);

        if (tolerancia(R, x)) {
            freeMatriz(&R);
            break;
        }

        freeMatriz(&R);
    }

    free(x);

    freeMatriz(&Jx);
    freeMatriz(&Fx);
}

void fP(double resultado[], rede r) {
    int iFP = 0;
    double angulo;

    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo == B_SWING)
            continue;

        resultado[iFP] = r.barras[j].tensao * r.barras[j].tensao * r.mNodal.condutancia.elemento[j][j];
        for (int k = 0; k < r.numBarras; k++) {
            if (k == j)
                continue;

            angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;

            resultado[iFP] += r.barras[j].tensao * r.barras[k].tensao * (r.mNodal.condutancia.elemento[j][k]*cos(angulo)
                                                                         - r.mNodal.susceptancia.elemento[j][k]*sin(angulo));
        }

        if (r.barras[j].tipo == B_PV)
            resultado[iFP] -= r.barras[j].potenciaAtiva;

        iFP++;
    }
}

void fQ(double resultado[], rede r) {
    int iFP = 0;
    double angulo;

    for (int j = 0; j < r.numBarras; j++) {
        if (r.barras[j].tipo != B_PQ)
            continue;

        resultado[iFP] = -r.barras[j].tensao * r.barras[j].tensao * r.mNodal.susceptancia.elemento[j][j];
        for (int k = 0; k < r.numBarras; k++) {
            if (k == j)
                continue;

            angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
            resultado[iFP] -= r.barras[j].tensao * r.barras[k].tensao * (r.mNodal.condutancia.elemento[j][k]*sin(angulo)
                                                                         + r.mNodal.susceptancia.elemento[j][k]*cos(angulo));
        }

        iFP++;
    }
}

void atualizaRede(double x[], rede r) {
    int jMatriz = 0, jPQ = 0;

    for (int k = 0; k < r.numBarras; k++) {
        if (r.barras[k].tipo == B_SWING)
            continue;

        r.barras[k].anguloTensao = x[jMatriz];

        if (r.barras[k].tipo == B_PQ) {
            r.barras[k].tensao = x[r.numPQ + r.numPV + jPQ];
            jPQ++;
        }

        jMatriz++;
    }
}

matriz funcaoDesvio(matriz M, rede r) {
    double *resultado = calloc(2*r.numPQ + r.numPV, sizeof(double));

    fP(resultado, r);
    fQ(&resultado[r.numPQ + r.numPV], r);

    for (int i = 0; i < 2*r.numPQ + r.numPV; i++)
        M.elemento[i][0] = -resultado[i];

    free(resultado);

    return M;
}


//double *theta = x;
//double *tensao = &x[r.numPQ + r.numPV];
matriz jacobianaDesvios(matriz M, rede r) {
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

            if (j == k) {
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

            if (j == k) {

                M.elemento[iFP][r.numPQ + r.numPV + jMatriz] = 2*r.barras[j].tensao * r.mNodal.condutancia.elemento[j][j];
                for (int i = 0; i < r.numBarras; i++) {
                    if (i == j)
                        continue;

                    angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                    M.elemento[iFP][r.numPQ + r.numPV + jMatriz] += r.barras[i].tensao * (r.mNodal.condutancia.elemento[j][i]*cos(angulo) - r.mNodal.susceptancia.elemento[j][i]*sin(angulo));
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

            if (j == k) {
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

            if (j == k) {
                M.elemento[r.numPQ + r.numPV + iFP][r.numPQ + r.numPV + jMatriz] = 0.0;

                for (int i = 0; i < r.numBarras; i++) {
                    if(i == j){
                        M.elemento[r.numPQ + r.numPV + iFP][r.numPQ + r.numPV + jMatriz] -= 2*r.barras[j].tensao * r.mNodal.susceptancia.elemento[j][j];
                    } else {
                        angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                        M.elemento[r.numPQ + r.numPV + iFP][r.numPQ + r.numPV + jMatriz] -= r.barras[i].tensao
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
