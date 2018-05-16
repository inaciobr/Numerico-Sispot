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
    int j = 0;
    for (int i = 0; i < r.numBarras; i++) {
        if (r.barras[i].tipo == B_SWING)
            continue;

        resultado[j] = r.barras[i].tensao;

        for (int k = 0; k < r.numBarras; k++) {
            double angulo = r.barras[k].anguloTensao - r.barras[i].anguloTensao;

            resultado[j] *= r.barras[k].tensao;
            resultado[j] *= r.mNodal.condutancia.elemento[i][k]*cos(angulo) - r.mNodal.susceptancia.elemento[i][k]*sin(angulo);
        }

        if (r.barras[i].tipo == B_PV)
            resultado[j] -= r.barras[i].potenciaAtiva;

        j++;
    }
}

void fQ(double resultado[], rede r) {
    int j = 0;
    for (int i = 0; i < r.numBarras; i++) {
        if (r.barras[i].tipo != B_PQ)
            continue;

        resultado[j] = -r.barras[i].tensao;

        for (int k = 0; k < r.numBarras; k++) {
            double angulo = r.barras[k].anguloTensao - r.barras[i].anguloTensao;

            resultado[j] *= r.barras[k].tensao;
            resultado[j] *= r.mNodal.condutancia.elemento[i][k] * sin(angulo) + r.mNodal.susceptancia.elemento[i][k] * cos(angulo);
        }

        j++;
    }
}

void funcaoDesvio(double x[], double resultado[], rede r) {
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

    fP(resultado, r);
    fQ(&resultado[r.numPQ + r.numPV], r);
}

matriz jacobianaDesvios(/*double x[], double res[], */rede r) {
    //double *theta = x;
    //double *tensao = &x[r.numPQ + r.numPV];
    matriz M = criaMatriz(2*r.numPQ + r.numPV, 2*r.numPQ + r.numPV);

    //del FP / del THETA => (r.numPQ + r.numPV) x (r.numPQ + r.numPV)

    int k,j = 0;
    double angulo;

     for (int l = 0; l < r.numBarras; l++) {
        if (r.barras[l].tipo == B_SWING)
            continue;

        k = 0;
        for (int i = 0; i < r.numBarras; i++) {
            if (r.barras[i].tipo == B_SWING)
                continue;

            if (j == k) {
                angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] +=r.barras[j].tensao*r.barras[i].tensao*(r.mNodal.condutancia.elemento[j][i]*sin(angulo)+r.mNodal.susceptancia.elemento[j][i]*cos(angulo));
            } else {
                angulo =r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] = -r.barras[j].tensao * r.barras[k].tensao*(r.mNodal.condutancia.elemento[j][k]*sin(angulo)+r.mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }

            k++;
        }

        j++;
    }


    //del FP / del V => (r.numPQ + r.numPV) x (r.numPQ)
    //j deve zerar

    int kPV = k;//ultimo valor da matriz dfp/dO + 1 na coluna
    int jQ = j;// ultimo valor da matriz dfp/dO + 1 na linha
    j = 0;

     for (int l = 0; l < r.numBarras; l++) {
        if (r.barras[l].tipo != B_PQ)
            continue;

        k = kPV;
        for (int i = 0; i < r.numBarras; i++) {
            if (r.barras[i].tipo != B_PQ)
                continue;

            if (j == k) {
                angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] += r.barras[j].tensao * r.barras[i].tensao*(r.mNodal.condutancia.elemento[j][i]*cos(angulo)-r.mNodal.susceptancia.elemento[j][i]*sin(angulo));
            } else {
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] = r.barras[j].tensao*(r.mNodal.condutancia.elemento[j][k]*cos(angulo)-r.mNodal.susceptancia.elemento[j][k]*sin(angulo));
            }

            k++;
        }

        j++;
    }


    //del FQ / del THETA => (r.numPQ) x (r.numPQ + r.numPV)
    j = jQ;
    k =0;
    for (int l = 0; l < r.numBarras; l++) {
        if (r.barras[l].tipo != B_PQ)
            continue;

        k = 0;
        for (int i = 0; i < r.numBarras; i++) {
            if (r.barras[i].tipo != B_PQ) //mudei de swing pra pq
                continue;

            if (j == k) {
                angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] +=r.barras[j].tensao*r.barras[i].tensao*(r.mNodal.condutancia.elemento[j][i]*cos(angulo)-r.mNodal.susceptancia.elemento[j][i]*sin(angulo));
            } else {
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] = -r.barras[j].tensao * r.barras[k].tensao;//*(r.mNodal.condutancia.elemento[j][k]*sin(angulo)+r.mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }//PROBLEMA//////////\//\/\\\/\//\/\/\

            k++;
        }

        j++;
    }


    //del FQ / del V => (r.numPQ) x (r.numPQ)
    j = jQ;
     for (int l = 0; l < r.numBarras; l++) {
        if (r.barras[l].tipo != B_PQ)
            continue;

        k = kPV;
        for (int i = 0; i < r.numBarras; i++) {
            if (r.barras[i].tipo != B_PQ)
                continue;

            if (j == k) {
                angulo = r.barras[i].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] += r.barras[j].tensao * r.barras[i].tensao ;//* (r.mNodal.condutancia.elemento[j][i]*sin(angulo)+r.mNodal.susceptancia.elemento[j][i]*cos(angulo));
            } else { //PROBLEMAAAAAAAA
                angulo = r.barras[k].anguloTensao - r.barras[j].anguloTensao;
                M.elemento[j][k] = -r.barras[j].tensao;//*(r.mNodal.condutancia.elemento[j][k]*sin(angulo)+r.mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }//PROBLEMAAAAA

            k++;
        }

        j++;
    }


    return M;

}

