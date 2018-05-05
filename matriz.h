#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int numLinhas;
    int numColunas;
    double **elemento;
} matriz;

//void criaMatriz();
double Det(matriz** matriz);
matriz** Inversa(matriz** matriz);
matriz** Cofatores(matriz** matriz);
void DecomposicaoLU(matriz** A,matriz** L,matriz** U);

