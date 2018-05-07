#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int numLinhas;
    int numColunas;
    double **elemento;
} matriz;

matriz criaMatriz(int numLinhas, int numColunas);
matriz copiaMatriz(matriz M);
void freeMatriz(matriz M);

matriz multiplicaConstante(matriz matriz, double constante);
matriz produtoMatriz(matriz M, matriz v);
matriz inversa(matriz M);
matriz transposta(matriz M);

double det(matriz A);

matriz decomposicaoLU(matriz A, int **permutacoes);
matriz permutaLinhasMatriz(matriz A, int *P);
matriz matrizCofatores(matriz A, int n, int m);
matriz cofatores(matriz M);
