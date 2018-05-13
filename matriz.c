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

void freeMatriz(matriz *A) {
    for (int i = 0; i < A->numLinhas; i++)
        free(A->elemento[i]);

    free(A->elemento);

    A->numLinhas = 0;
    A->numColunas = 0;
}

matriz multiplicaConstante(matriz A, double constante) {
    for (int i = 0; i < A.numLinhas; i++)
        for (int j = 0; j < A.numColunas; j++)
            A.elemento[i][j] *= constante;

    return A;
}

double det(matriz M) {
    if (M.numColunas != M.numLinhas)
        return 0;

    int *P;
    matriz LU = decomposicaoLU(M, &P);

    int N = P[M.numLinhas];
    free(P);

    double det = 1;
    for (int i = 0; i < M.numLinhas; i++)
        det *= LU.elemento[i][i];

    freeMatriz(&LU);

    return det * (N%2 ? -1 : 1);
}

matriz transposta(matriz M) {
    matriz matrizT = criaMatriz(M.numColunas, M.numLinhas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
                matrizT.elemento[j][i] = M.elemento[i][j];

    freeMatriz(&M);

    return matrizT;
}

matriz copiaMatriz(matriz M) {
    matriz copiaM = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            copiaM.elemento[i][j] = M.elemento[i][j];

    return copiaM;
}

matriz decomposicaoLU(matriz M, int **permutacoes) {
    matriz LU = copiaMatriz(M);
    int *P = malloc((M.numLinhas + 1) * sizeof(double));
    P[M.numLinhas] = 0;

    for (int k = 0; k < M.numLinhas; k++) {
        double maxPivo = -1.0;

        for (int i = k; i < M.numLinhas; i++) {
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

            P[M.numLinhas]++;
        }

        for (int j = k + 1; j < M.numLinhas; j++) {
            for (int i = 0; i <= k - 1; i++)
                LU.elemento[k][j] -= LU.elemento[k][i] * LU.elemento[i][j];

            LU.elemento[j][k] /= LU.elemento[k][k];
        }
    }

    *permutacoes = P;
    return LU;
}

matriz permutaLinhasMatriz(matriz A, int P[]) {
    for (int i = 0; i < A.numLinhas; i++)
        if (i != P[i]) {
            double *temp = A.elemento[i];
            A.elemento[i] = A.elemento[P[i]];
            A.elemento[P[i]] = temp;
        }

    return A;
}

double cofator(matriz M, int linha, int coluna) {
    matriz mCofatores = criaMatriz(M.numLinhas - 1, M.numColunas - 1);
    double resultado;

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            if (i != linha && j != coluna)
                mCofatores.elemento[(i <= linha ? i : i - 1)][(j <= coluna ? j : j - 1)] = M.elemento[i][j];

    resultado = ((linha + coluna)%2 ? -1 : 1) * det(mCofatores);

    freeMatriz(&mCofatores);

    return resultado;
}

matriz matrizCofatores(matriz M) {
    matriz cofatores = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            cofatores.elemento[i][j] = cofator(M, i, j);

    return cofatores;
}

matriz inversa(matriz M) {
    if (M.numColunas != M.numLinhas)
        return M;

    matriz cof = matrizCofatores(M);
    matriz inversa = multiplicaConstante(transposta(cof), 1/det(M));

    return inversa;
}

matriz produtoMatriz(matriz M1, matriz M2){
    matriz produto = criaMatriz(M1.numLinhas, M2.numColunas);

    for (int i = 0; i < M1.numLinhas; i++)
        for (int j = 0; j < M2.numColunas; j++)
                for (int k = 0; k < M1.numColunas; k++)
                    produto.elemento[i][j] += M1.elemento[i][k] * M2.elemento[k][j];

    return produto;
}

/// MUDAR PARA LU
/// Talvez resultado em double (?)
matriz resolveSistemaLinear(matriz A, matriz b) {
    matriz x = produtoMatriz(inversa(A), b);

    return x;
}

void printMatriz(matriz M) {
    for (int i = 0; i < M.numLinhas; i++) {
        for (int j = 0; j < M.numColunas; j++)
            printf("%15.8f ", M.elemento[i][j]);

        printf("\n");
    }
}

