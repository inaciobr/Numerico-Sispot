#include "matriz.h"

matriz criaMatriz(int numLinhas, int numColunas) {
    matriz M;
    M.numLinhas = numLinhas;
    M.numColunas = numColunas;

    M.elemento = malloc(numLinhas * sizeof(double *));

    for (int i = 0; i < numLinhas; i++)
        M.elemento[i] = calloc(numColunas, sizeof(double));

    return M;
}

void freeMatriz(matriz M) {
    for (int i = 0; i < M.numLinhas; i++) {
        free(M.elemento[i]);
    }

    free(M.elemento);
}

matriz multiplicaConstante(matriz M, double constante) {
    matriz resultado = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            resultado.elemento[i][j] = constante * M.elemento[i][j];

    return resultado;
}

double det(matriz A) {
    int *P;
    matriz LU = decomposicaoLU(A, &P);
    int N = P[A.numLinhas];

    double det = 1;
    for (int i = 0; i < A.numLinhas; i++)
        det *= LU.elemento[i][i];

    freeMatriz(LU);
    free(P);

    return det * (N%2 ? -1 : 1);
}

matriz transposta(matriz M) {
    matriz matrizT = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
                matrizT.elemento[i][j] = M.elemento[j][i];

    return matrizT;
}

matriz copiaMatriz(matriz M) {
    matriz copiaM = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            copiaM.elemento[i][j] = M.elemento[i][j];

    return copiaM;
}

matriz decomposicaoLU(matriz A, int **permutacoes) {
    matriz LU = copiaMatriz(A);
    int *P = malloc((A.numLinhas + 1) * sizeof(double));
    P[A.numLinhas] = 0;

    for (int k = 0; k < A.numLinhas; k++) {
        double maxPivo = -1.0;

        for (int i = k; i < A.numLinhas; i++) {
            for (int j = 0; j <= k - 1; j++)
                LU.elemento[i][k] -= LU.elemento[i][j] * LU.elemento[j][k];

            if (abs(LU.elemento[i][k]) > maxPivo) {
                maxPivo = LU.elemento[i][k];
                P[k] = i;
            }
        }

        if (P[k] != k) {
            double *temp = LU.elemento[k];
            LU.elemento[k] = LU.elemento[P[k]];
            LU.elemento[P[k]] = temp;

            P[A.numLinhas]++;
        }

        for (int j = k + 1; j < A.numLinhas; j++) {
            for (int i = 0; i <= k - 1; i++)
                LU.elemento[k][j] -= LU.elemento[k][i] * LU.elemento[i][j];

            LU.elemento[j][k] /= LU.elemento[k][k];
        }
    }

    *permutacoes = P;
    return LU;
}

matriz permutaLinhasMatriz(matriz A, int *P) {
    matriz M = copiaMatriz(A);

    for (int i = 0; i < M.numLinhas; i++) {
        if (i != P[i]) {
            double *temp = M.elemento[i];
            M.elemento[i] = M.elemento[P[i]];
            M.elemento[P[i]] = temp;
        }
    }

    return M;
}

matriz matrizCofatores(matriz M, int linha, int coluna) {
    matriz MCofatores = criaMatriz(M.numLinhas - 1, M.numColunas - 1);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            if (i != linha && j != coluna)
                MCofatores.elemento[(i <= linha ? i : i - 1)][(j <= coluna ? j : j - 1)] = M.elemento[i][j];

    return MCofatores;
}

matriz cofatores(matriz M) {
    matriz cofatores = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++) {
        for (int j = 0; j < M.numColunas; j++) {
            matriz cofij = matrizCofatores(M, i, j);
            cofatores.elemento[i][j] = ((i + j)%2 ? -1 : 1) * det(cofij);
            freeMatriz(cofij);
        }
    }

    return cofatores;
}

matriz inversa(matriz M) {
    matriz cof = cofatores(M);
    matriz inversa = multiplicaConstante(transposta(cofatores(M)), 1/det(M));
    freeMatriz(cof);

    return inversa;
}

matriz produtoMatriz(matriz A, matriz B){
    matriz produto = criaMatriz(A.numLinhas, B.numColunas);

    for (int i = 0; i < A.numLinhas ; i++)
        for (int j = 0; j < B.numColunas ; j++)
                for (int k = 0; k < A.numColunas; k++)
                    produto.elemento[i][j] += A.elemento[i][k] * B.elemento[k][j];

    return produto;
}
