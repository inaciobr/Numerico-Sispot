#ifndef MATRIZ_H_INCLUDED
#define MATRIZ_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef struct {
    int numLinhas;      /** Número de linhas da matriz. */
    int numColunas;     /** Número de colunas da matriz. */
    double **elemento;  /** Matriz do tipo double onde elementos da matriz são armazenados. */
} matriz;

matriz criaMatriz(int numLinhas, int numColunas);
matriz copiaMatriz(matriz M);
void freeMatriz(matriz *M);

matriz multiplicaConstante(matriz A, double constante);
matriz produtoMatriz(matriz M, matriz v);
matriz inversa(matriz M);
matriz transposta(matriz M);

double det(matriz A);
double cofator(matriz A, int n, int m);

matriz decomposicaoLU(matriz A, int **permutacoes);
matriz resolveSistemaLinear(matriz A, matriz b);
matriz cofatores(matriz M);

void permutaLinhasMatriz(matriz A, int P[]);
void printMatriz(matriz M);

#endif // MATRIZ_H_INCLUDED
