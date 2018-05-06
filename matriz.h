#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int numLinhas;
    int numColunas;
    double **elemento;
} matriz;

matriz* criaMatriz();
matriz *criaMatrizCof(matriz A, int n, int m);
matriz* multiplicaConstante(double cnst, matriz matriz);
void DecomposicaoLU(matriz A,matriz* L,matriz* U);
double Det(matriz matriz);
matriz* Cofatores(matriz matriz);
matriz* Transposta(matriz matriz);
matriz* Inversa(matriz matriz);
