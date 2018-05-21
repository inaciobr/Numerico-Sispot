#include "rede.h"

/**
 *
 */
rede* leituraRede(char arquivo[]) {
    rede *r = malloc(sizeof(rede));

    r->numPQ = 0;
    r->numPV = 0;
    r->numSwing = 0;

    memcpy(r->nome, arquivo, strlen(arquivo) + 1);

	r->potenciaAtivaAbsorvida = 0.0;
	r->potenciaAtivaGerada = 0.0;
	r->perdaAtiva = 0.0;

    char caminho1[256], caminho2[256];

    sprintf(caminho1, "%s_DadosBarras.txt", arquivo);
    sprintf(caminho2, "%s_YNodal.txt", arquivo);

    leituraBarra(r, caminho1);

	if (r->barras == NULL) {
		free(r);
		return NULL;
	}

	leituraNodal(r, caminho2);

    return r;
}

/**
 *
 */
void leituraBarra(rede *r, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("N�o foi poss�vel ler o arquivo %s.", arquivo);
		r->barras = NULL;
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
			r->numPQ++;
            fscanf(fp, "%lf %lf", &barrasArquivo[i].potenciaAtivaNominal,
                                  &barrasArquivo[i].potenciaReativaNominal);

            barrasArquivo[i].anguloTensao = 0.0;
            barrasArquivo[i].tensao = barrasArquivo[i].tensaoNominal;

            break;

        case B_PV:
			r->numPV++;
            fscanf(fp, "%lf %lf", &barrasArquivo[i].potenciaAtivaNominal,
                                  &barrasArquivo[i].tensao);

            barrasArquivo[i].anguloTensao = 0.0;

            break;

        case B_SWING:
			r->numSwing++;
            fscanf(fp, "%lf %lf", &barrasArquivo[i].tensao,
                                  &barrasArquivo[i].anguloTensao);

            barrasArquivo[i].potenciaAtivaNominal = 0.0;
            barrasArquivo[i].potenciaReativaNominal = 0.0;

            break;
        }
    }

    fclose(fp);
    r->barras = barrasArquivo;
}

/**
 *
 */
void leituraNodal(rede *r, char arquivo[]) {
    FILE *fp;
    fp = fopen(arquivo, "r");
    if (fp == NULL) {
        printf("N�o foi poss�vel ler o arquivo %s.", arquivo);
        return;
    }

    int N;
    fscanf(fp, "%d", &N);

    r->mNodal.condutancia = criaMatriz(r->numBarras, r->numBarras);
    r->mNodal.susceptancia = criaMatriz(r->numBarras, r->numBarras);

    int linha, coluna;
    for (int i = 0; i < N; i++) {
        fscanf(fp, "%d %d", &linha, &coluna);
        fscanf(fp, "%lf %lf", &r->mNodal.condutancia.elemento[linha][coluna],
                              &r->mNodal.susceptancia.elemento[linha][coluna]);
    }

    fclose(fp);
}

/**
 *
 */
void freeRede(rede *r) {
    r->potenciaAtivaGerada = r->potenciaAtivaAbsorvida = r->perdaAtiva = 0.0;
    r->numBarras = r->numPQ = r->numPV = r->numSwing = 0;

    freeMatriz(&r->mNodal.condutancia);
    freeMatriz(&r->mNodal.susceptancia);
    free(r->barras);
}

/**
 *
 */
void fluxoDePotenciaNewton(rede *r) {
    matriz Jx = criaMatriz(2*r->numPQ + r->numPV, 2*r->numPQ + r->numPV);
    matriz Fx = criaMatriz(2*r->numPQ + r->numPV, 1);
    matriz R;

    double *x = calloc(2*r->numPQ + r->numPV, sizeof(double));

    // Valores iniciais de x.
    int jMatriz = 0, jPQ = 0;
    for (int k = 0; k < r->numBarras; k++) {
        if (r->barras[k].tipo == B_SWING)
            continue;

        x[jMatriz] = r->barras[k].anguloTensao;

        if (r->barras[k].tipo == B_PQ) {
            x[r->numPQ + r->numPV + jPQ] = r->barras[k].tensao;
            jPQ++;
        }

        jMatriz++;
    }
	int k;
    for (k = 0; k < MAX_ITERACOES; k++) {
        Fx = funcaoDesvio(Fx, r);
        Jx = jacobianaDesvios(Jx, r);

        R = resolveSistemaLinear(Jx, Fx);

        for (int i = 0; i < R.numLinhas; i++)
            x[i] += R.elemento[i][0];

		atualizaBarrasX(x, r);

        if (tolerancia(R, x)) {
            freeMatriz(&R);
            break;
        }

        freeMatriz(&R);
    }

    free(x);

    freeMatriz(&Jx);
    freeMatriz(&Fx);

	atualizaRedePU(r);
}

/**
 *
 */
void fP(double resultado[], rede *r) {
    int iFP = 0;
    double angulo;

    for (int j = 0; j < r->numBarras; j++) {
        if (r->barras[j].tipo == B_SWING)
            continue;

        resultado[iFP] = r->barras[j].tensao * r->barras[j].tensao * r->mNodal.condutancia.elemento[j][j];
        for (int k = 0; k < r->numBarras; k++) {
            if (k == j)
                continue;

            angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;

            resultado[iFP] += r->barras[j].tensao * r->barras[k].tensao * (r->mNodal.condutancia.elemento[j][k]*cos(angulo)
                                                                         - r->mNodal.susceptancia.elemento[j][k]*sin(angulo));
        }

        if (r->barras[j].tipo == B_PV)
            resultado[iFP] -= r->barras[j].potenciaAtivaNominal;

        iFP++;
    }
}

/**
 *
 */
void fQ(double resultado[], rede *r) {
    int iFP = 0;
    double angulo;

    for (int j = 0; j < r->numBarras; j++) {
        if (r->barras[j].tipo != B_PQ)
            continue;

        resultado[iFP] = -r->barras[j].tensao * r->barras[j].tensao * r->mNodal.susceptancia.elemento[j][j];
        for (int k = 0; k < r->numBarras; k++) {
            if (k == j)
                continue;

            angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;
            resultado[iFP] -= r->barras[j].tensao * r->barras[k].tensao * (r->mNodal.condutancia.elemento[j][k]*sin(angulo)
                                                                         + r->mNodal.susceptancia.elemento[j][k]*cos(angulo));
        }

        iFP++;
    }
}

/**
 *
 */
void atualizaBarrasX(double x[], rede *r) {
    int jMatriz = 0, jPQ = 0;

    for (int k = 0; k < r->numBarras; k++) {
        if (r->barras[k].tipo == B_SWING)
            continue;

        r->barras[k].anguloTensao = x[jMatriz];

        if (r->barras[k].tipo == B_PQ) {
            r->barras[k].tensao = x[r->numPQ + r->numPV + jPQ];
            jPQ++;
        }

        jMatriz++;
    }
}

/**
 *
 */
void atualizaRedePU(rede *r) {
	r->perdaAtiva = r->potenciaAtivaGerada = r->potenciaAtivaAbsorvida = 0.0;
	double angulo, condutanciaCarga;

	for (int j = 0; j < r->numBarras; j++) {
		r->barras[j].valorPorUnidade = r->barras[j].tensao / r->barras[j].tensaoNominal;

		condutanciaCarga = 0.0;
		for (int k = 0; k < r->numBarras; k++) {
			angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;
			r->potenciaAtivaGerada += 3 * r->barras[j].tensao * r->barras[k].tensao * (r->mNodal.condutancia.elemento[j][k]*cos(angulo) - r->mNodal.susceptancia.elemento[j][k]*sin(angulo));

            condutanciaCarga += r->mNodal.condutancia.elemento[j][k];
		}

        r->barras[j].potenciaAtiva = 0.0;
        r->barras[j].potenciaReativa = 0.0;

		r->potenciaAtivaAbsorvida += 3 * condutanciaCarga * r->barras[j].tensao * r->barras[j].tensao;
	}

	r->perdaAtiva = r->potenciaAtivaGerada - r->potenciaAtivaAbsorvida;
}

/**
 *
 */
matriz funcaoDesvio(matriz M, rede *r) {
    double *resultado = calloc(2*r->numPQ + r->numPV, sizeof(double));

    fP(resultado, r);
    fQ(&resultado[r->numPQ + r->numPV], r);

    for (int i = 0; i < 2*r->numPQ + r->numPV; i++)
        M.elemento[i][0] = -resultado[i];

    free(resultado);

    return M;
}


//double *theta = x;
//double *tensao = &x[r->numPQ + r->numPV];
/**
 *
 */
matriz jacobianaDesvios(matriz M, rede *r) {
    double angulo;

    //del FP / del THETA => (r->numPQ + r->numPV) x (r->numPQ + r->numPV)
    int jMatriz = 0, iFP = 0;
    for (int j = 0; j < r->numBarras; j++) {
        if (r->barras[j].tipo == B_SWING)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r->numBarras; k++) {
            if (r->barras[k].tipo == B_SWING)
                continue;

            if (j == k) {
                M.elemento[iFP][jMatriz] = 0.0;

                for (int i = 0; i < r->numBarras; i++) {
                    if(i == j)
                        continue;
                    angulo = r->barras[i].anguloTensao - r->barras[j].anguloTensao;
                    M.elemento[iFP][jMatriz] += r->barras[i].tensao * (r->mNodal.condutancia.elemento[j][i]*sin(angulo) + r->mNodal.susceptancia.elemento[j][i]*cos(angulo));
                }

                M.elemento[iFP][jMatriz] *= r->barras[j].tensao;
            } else {
                angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;
                M.elemento[iFP][jMatriz] = -r->barras[j].tensao * r->barras[k].tensao * (r->mNodal.condutancia.elemento[j][k]*sin(angulo) + r->mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }

    //del FP / del V => (r->numPQ + r->numPV) x (r->numPQ)
    iFP = 0;
    for (int j = 0; j < r->numBarras; j++) {
        if (r->barras[j].tipo == B_SWING)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r->numBarras; k++) {
            if (r->barras[k].tipo != B_PQ)
                continue;

            if (j == k) {

                M.elemento[iFP][r->numPQ + r->numPV + jMatriz] = 2*r->barras[j].tensao * r->mNodal.condutancia.elemento[j][j];
                for (int i = 0; i < r->numBarras; i++) {
                    if (i == j)
                        continue;

                    angulo = r->barras[i].anguloTensao - r->barras[j].anguloTensao;
                    M.elemento[iFP][r->numPQ + r->numPV + jMatriz] += r->barras[i].tensao * (r->mNodal.condutancia.elemento[j][i]*cos(angulo) - r->mNodal.susceptancia.elemento[j][i]*sin(angulo));
                }

            } else {
                angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;
                M.elemento[iFP][r->numPQ + r->numPV + jMatriz] = r->barras[j].tensao * (r->mNodal.condutancia.elemento[j][k]*cos(angulo) - r->mNodal.susceptancia.elemento[j][k]*sin(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }


    //del FQ / del THETA => (r->numPQ) x (r->numPQ + r->numPV)
    iFP = 0;
    for (int j = 0; j < r->numBarras; j++) {
        if (r->barras[j].tipo != B_PQ)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r->numBarras; k++) {
            if (r->barras[k].tipo == B_SWING)
                continue;

            if (j == k) {
                M.elemento[r->numPQ + r->numPV + iFP][jMatriz] = 0.0;

                for (int i = 0; i < r->numBarras; i++) {
                    if(i == j)
                        continue;

                    angulo = r->barras[i].anguloTensao - r->barras[j].anguloTensao;
                    M.elemento[r->numPQ + r->numPV + iFP][jMatriz] += r->barras[i].tensao*(r->mNodal.condutancia.elemento[j][i]*cos(angulo) - r->mNodal.susceptancia.elemento[j][i]*sin(angulo));
                }
                M.elemento[r->numPQ + r->numPV + iFP][jMatriz] *= r->barras[j].tensao;
            } else {
                angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;
                M.elemento[r->numPQ + r->numPV + iFP][jMatriz] = -r->barras[j].tensao * r->barras[k].tensao * (r->mNodal.condutancia.elemento[j][k]*cos(angulo) - r->mNodal.susceptancia.elemento[j][k]*sin(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }


    //del FQ / del V => (r->numPQ) x (r->numPQ)
    iFP = 0;
    for (int j = 0; j < r->numBarras; j++) {
        if (r->barras[j].tipo != B_PQ)
            continue;

        jMatriz = 0;
        for (int k = 0; k < r->numBarras; k++) {
            if (r->barras[k].tipo != B_PQ)
                continue;

            if (j == k) {
                M.elemento[r->numPQ + r->numPV + iFP][r->numPQ + r->numPV + jMatriz] = 0.0;

                for (int i = 0; i < r->numBarras; i++) {
                    if(i == j){
                        M.elemento[r->numPQ + r->numPV + iFP][r->numPQ + r->numPV + jMatriz] -= 2*r->barras[j].tensao * r->mNodal.susceptancia.elemento[j][j];
                    } else {
                        angulo = r->barras[i].anguloTensao - r->barras[j].anguloTensao;
                        M.elemento[r->numPQ + r->numPV + iFP][r->numPQ + r->numPV + jMatriz] -= r->barras[i].tensao * (r->mNodal.condutancia.elemento[j][i]*sin(angulo) + r->mNodal.susceptancia.elemento[j][i]*cos(angulo));
                    }
                }
            } else {
                angulo = r->barras[k].anguloTensao - r->barras[j].anguloTensao;
                M.elemento[r->numPQ + r->numPV + iFP][r->numPQ + r->numPV + jMatriz] = -r->barras[j].tensao * (r->mNodal.condutancia.elemento[j][k]*sin(angulo) + r->mNodal.susceptancia.elemento[j][k]*cos(angulo));
            }

            jMatriz++;
        }

        iFP++;
    }

    return M;
}

double perdaTrecho(rede *r, int barra1, int barra2) {
    double complex tensao1 = r->barras[barra1].tensao * cexp(1i*r->barras[barra1].anguloTensao);
    double complex tensao2 = r->barras[barra2].tensao * cexp(1i*r->barras[barra2].anguloTensao);
    double deltaV = cabs(tensao1 - tensao2);

    return -3 * deltaV * deltaV * r->mNodal.condutancia.elemento[barra1][barra2];
}

double fluxoPotencia(rede *r, int barra1, int barra2) {
    double complex tensao1 = r->barras[barra1].tensao * cexp(1i*r->barras[barra1].anguloTensao);
    double complex tensao2 = r->barras[barra2].tensao * cexp(1i*r->barras[barra2].anguloTensao);
    double complex admitancia = -r->mNodal.condutancia.elemento[barra1][barra2] - 1i*r->mNodal.susceptancia.elemento[barra1][barra2];

    return 3 * creal(tensao1 * conj((tensao1 - tensao2) * admitancia));
}

/**
 *
 */
void printRede(rede *r, FILE *saida) {
	printf("Analise de: %s\n\n", r->nome);

	fprintf(saida, "                       Resultados das barras\n");
	fprintf(saida, "Barra | Modulo (PU) |  Angulo (o) | Modulo da tensao (V)\n");

	for (int k = 0; k < r->numBarras; k++)
		fprintf(saida, "%5d | %11.6f | %11.4f | %15.3f\n", r->barras[k].id, r->barras[k].valorPorUnidade, rad2Graus(r->barras[k].anguloTensao), r->barras[k].tensao);

    fprintf(saida, "\n\n");

	fprintf(saida, "                       Resultados dos trechos\n");
	fprintf(saida, "Barra inicial | Barra final | Potencia ativa (kW) | Perda ativa (kW)\n");

    for (int k = 0; k < r->numBarras; k++)
        for (int j = 0; j < r->numBarras; j++)
            if ((r->mNodal.condutancia.elemento[k][j] || r->mNodal.susceptancia.elemento[k][j]) && k != j)
                fprintf(saida, "%13d | %11d | %19.3f | %16.3f\n", r->barras[k].id,
                                                                           r->barras[j].id,
                                                                           fluxoPotencia(r, k, j) / 1000.,
                                                                           perdaTrecho(r, k, j) / 1000.);

    fprintf(saida, "\n\n");

	fprintf(saida, "                       Resultados globais\n");
	fprintf(saida, "Potencia ativa total gerada:   %15.3f (kW)\n", r->potenciaAtivaGerada / 1000.);
	fprintf(saida, "Potencia ativa total de carga: %15.3f (kW)\n", r->potenciaAtivaAbsorvida / 1000.);
	fprintf(saida, "Perda ativa total:             %15.3f (kW)\n", r->perdaAtiva / 1000.);

}

/**
 *
 */
void printDadosRede(rede *r) {
    printRede(r, stdout);
}

/**
 *
 */
void arquivarDadosRede(rede *r) {
    FILE *fp;
    char caminho[256];

    sprintf(caminho, "%s_resultados.txt", r->nome);
    fp = fopen(caminho, "w");

    printRede(r, fp);
    printf("Dados da rede salvos em %s_resultados.txt\n", r->nome);

    fclose(fp);
}
