#ifndef BARRA_H_INCLUDED
#define BARRA_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

enum tipoBarra {
    B_PQ,         // Tipo 0
    B_PV,         // Tipo 1
    B_SWING       // Tipo 2
};

typedef struct {
    int id;
    int tipo;
    double tensao;
    double anguloTensao;
    double potenciaAtiva;
    double potenciaReativa;
} barra;

int leituraBarra(barra** barras, char arquivo[]) ;

#endif // BARRA_H_INCLUDED
