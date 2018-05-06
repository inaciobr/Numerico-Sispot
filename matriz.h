#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int numLinhas;
    int numColunas;
    double **elemento;
} matriz;

matriz criaMatriz(int numLinhas, int numColunas);
void freeMatriz(matriz *M);
double det(matriz A);


void decomposicaoLU(matriz A, matriz *L, matriz *U);
matriz *criaMatrizCof(matriz A, int n, int m);
matriz* multiplicaConstante(double cnst, matriz matriz);
matriz* Cofatores(matriz matriz);
matriz* Transposta(matriz matriz);
matriz* Inversa(matriz matriz);
