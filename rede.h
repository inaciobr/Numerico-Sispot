#ifndef REDE_H_INCLUDED
#define REDE_H_INCLUDED

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matriz.h"
#include "numerico.h"

typedef struct {
    matriz condutancia;
    matriz susceptancia;
} nodal;

typedef struct {
    int id;
    int tipo;

    double tensaoNominal;
	double potenciaAtivaNominal;
	double potenciaReativaNominal;

    double valorPorUnidade;
    double tensao;
	double anguloTensao;
} barra;

typedef struct {
    char nome[256];

	double potenciaAtivaGerada;
	double potenciaAtivaAbsorvida;
	double perdaAtiva;

    nodal mNodal;

    int numBarras;
    barra *barras;

    int numPQ;
    int numPV;
    int numSwing;
} rede;

enum tipoBarra {
    B_PQ,         // Tipo 0
    B_PV,         // Tipo 1
    B_SWING       // Tipo 2
};

rede* leituraRede(char arquivo[]);
void leituraBarra(rede *r, char arquivo[]);
void leituraNodal(rede *r, char arquivo[]);

void freeRede(rede *r);

int fluxoDePotenciaNewton(rede *r);

matriz funcaoDesvio(matriz M, rede *r);
matriz jacobianaDesvios(matriz M, rede *r);

void fP(double resultado[], rede *r);
void fQ(double resultado[], rede *r);

void atualizaBarrasX(double x[], rede *r);
void atualizaRede(rede *r);

double perdaTrecho(rede *r, int barra1, int barra2);

void printDadosRede(rede *r);
void arquivarDadosRede(rede *r);

#endif // REDE_H_INCLUDED
