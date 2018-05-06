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
            resultado.elemento[i][j] = constante * (M.elemento[i][j]);

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

matriz *criaMatrizCof(matriz M, int n, int m) {
    matriz *Mcof;
    double** elemento = malloc((M.numColunas-1)* sizeof(double*));

    for(int i = 0;i < M.numColunas - 1; i++)
        elemento[i] = malloc((M.numLinhas-1) * sizeof(double));

    for(int i = 0;i < n;i++){
        for(int j = 0;j < m ;j++){
            elemento[i][j] = M.elemento[i][j];
        }
    }

    for(int i = n+1;i < M.numLinhas ;i++){
        for(int j = 0;j < m ;j++){
            elemento[i-1][j] = M.elemento[i][j];
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = m;j < M.numColunas;j++){
            elemento[i][j-1] = M.elemento[i][j];
        }
    }

    for(int i = n;i < M.numLinhas;i++){
        for(int j = m;j < M.numColunas;j++){
            elemento[i-1][j-1] = M.elemento[i][j];
        }
    }

    Mcof->elemento = elemento;
    return Mcof;
}

/****TA DANDO MELDA***/
matriz* Cofatores(matriz M){
    matriz* cof, *cofij;
    double** elemento = malloc(M.numColunas* sizeof(double*));

    for(int i = 0;i < M.numColunas; i++)
        elemento[i] = malloc(M.numLinhas*sizeof(double));


    for(int i = 0;i < M.numLinhas;i++){
        for(int j = 0;j < M.numColunas;j++){
            cofij = criaMatrizCof(M,i,j);
            elemento[i][j] = pow((-1),(i+j))*det(*cofij);
        }
    }
    cof->elemento = elemento;
    return cof;
}

matriz* Inversa(matriz M){
    //return multiplicaConstante(det(M),*Transposta(*(Cofatores(M))));
}
