#ifndef REDE_H_INCLUDED
#define REDE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
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
    double anguloTensao;
    double potenciaAtiva;
    double potenciaReativa;

    double valorPorUnidade;
    double tensao;
} barra;

typedef struct {
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


rede leituraRede(char rede[]);
void leituraBarra(rede *r, char arquivo[]);
void leituraNodal(rede *r, char arquivo[]);

void freeNodal(nodal *mNodal);

void fluxoDePotenciaNewton(rede r);

matriz funcaoDesvio(matriz M, rede r);
matriz jacobianaDesvios(matriz M, rede r);

void fP(double resultado[], rede r);
void fQ(double resultado[], rede r);

void atualizaRede(double x[], rede r);

#endif // REDE_H_INCLUDED