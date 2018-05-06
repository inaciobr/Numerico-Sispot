#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int numLinhas;
    int numColunas;
    double **elemento;
} matriz;

matriz criaMatriz(int numLinhas, int numColunas);
void freeMatriz(matriz *M);
double det(matriz A);


void decomposicaoLU(matriz A, matriz *L, matriz *U);
matriz **Inversa(matriz **matriz);
matriz **Cofatores(matriz **matriz);
