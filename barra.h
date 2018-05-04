#include <stdio.h>

enum tipoBarra {
    PQ,         // Tipo 0
    PV,         // Tipo 1
    SWING       // Tipo 2
};

typedef struct {
    int id;
    int tipo;
    double tensaoNominal;
    double anguloTensao;
    double potenciaAtiva;
    double potenciaReativa;
} barra;

int leituraBarra(barra** barras, char arquivo[]) ;
