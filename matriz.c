#include "matriz.h"

matriz criaMatriz(int numLinhas, int numColunas) {
    matriz M;
    M.numLinhas = numLinhas;
    M.numColunas = numColunas;

    M.elemento = malloc(numLinhas * sizeof(double *));

    for (int i = 0; i < numLinhas; i++)
        M.elemento[i] = malloc(numColunas * sizeof(double));

    return M;
}

void freeMatriz(matriz *M) {
    for (int i = 0; i < M->numLinhas; i++) {
        free(M->elemento[i]);
    }

    free(M->elemento);
}

matriz multiplicaConstante(matriz M, double constante) {
    matriz resultado = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            resultado.elemento[i][j] = constante * M.elemento[i][j];

    return resultado;
}

double det(matriz A) {
    matriz L, U;
    decomposicaoLU(A, &L, &U);

    double det = 1;
    for (int i = 0; i < U.numColunas ;i++)
        det *= U.elemento[i][i];

    freeMatriz(&L);
    freeMatriz(&U);

    return det;
}

matriz transposta(matriz M){
    matriz matrizT = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
                matrizT.elemento[i][j] = M.elemento[j][i];

    return matrizT;
}

void decomposicaoLU(matriz A, matriz *L, matriz *U){
    *L = criaMatriz(A.numLinhas, A.numColunas);
    *U = criaMatriz(A.numLinhas, A.numColunas);

    for (int i = 0; i < A.numLinhas; i++) {
        for (int j = 0; j < A.numColunas; j++) {
            U->elemento[i][j] = A.elemento[i][j];
            if (i == j)
                L->elemento[i][j] = 1;
            else
                L->elemento[i][j] = 0;
        }
    }


    double *P = malloc(A.numLinhas * sizeof(double));
    for (int i = 0; i < A.numLinhas; i++)
        P[i] = i;

    for (int i = 0; i < A.numLinhas; i++) {
        double maxA = 0.0;
        int linhaMaxA = 0;

        for (int j = 0; j < U->numLinhas; j++) {
            // Encontra o maior pivô.
            if (abs(U->elemento[j][i]) > maxA) {
                maxA = abs(U->elemento[i][j]);
                linhaMaxA = j;
            }
        }

        if (maxA == 0)
            return;

        if (linhaMaxA != i) {
            // Troca linhas da matriz de permutação.
            double temp = P[i];
            P[i] = P[linhaMaxA];
            P[linhaMaxA] = temp;

            // Troca linhas da matriz de A.
            double *tempPtr = U->elemento[i];
            U->elemento[i] = U->elemento[linhaMaxA];
            U->elemento[linhaMaxA] = tempPtr;
        }

        for (int j = i + 1; j < A.numLinhas; j++) {
            L->elemento[j][i] = U->elemento[j][i] / U->elemento[i][i];

            for (int k = 1; k < A.numLinhas; k++) {
                U->elemento[j][k] -= L->elemento[j][i] * U->elemento[i][k];
            }
        }
    }

    for (int j = 0; j < A.numLinhas; j++) {
        for (int k = 0; k < j; k++) {
            U->elemento[j][k] = 0;
        }
    }
}

matriz matrizCofatores(matriz M, int linha, int coluna) {
    matriz MCofatores = criaMatriz(M.numLinhas - 1, M.numColunas - 1);

    for (int i = 0; i < M.numLinhas; i++)
        for(int j = 0; j < M.numColunas; j++)
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
            freeMatriz(&cofij);
        }
    }

    return cofatores;
}

matriz inversa(matriz M) {
    matriz cof = cofatores(M);
    matriz inversa = multiplicaConstante(transposta(cofatores(M)), 1/det(M));
    freeMatriz(&cof);

    return inversa;
}
