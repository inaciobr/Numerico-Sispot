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

matriz multiplicaConstante(matriz matriz, double constante);
matriz transposta(matriz M);

double det(matriz A);

void decomposicaoLU(matriz A, matriz *L, matriz *U);

matriz matrizCofatores(matriz A, int n, int m);
matriz cofatores(matriz M);
matriz inversa(matriz M);
