#include "matriz.h"

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
///TA TOP
matriz* multiplicaConstante(double cnst, matriz M){
    matriz *mult;
    double** elemento = malloc(M.numColunas* sizeof(double*));

    for(int i = 0;i < M.numColunas; i++)
        elemento[i] = malloc(M.numLinhas*sizeof(double));

    for(int i=0;i<M.numLinhas;i++)
        for(int j=0;j<M.numColunas;j++)
            elemento[i][j] = cnst * (M.elemento[i][j]);

    mult->elemento = elemento;
    return mult;
}

void DecomposicaoLU(matriz A,matriz* L,matriz* U){

}
///TA BIXANDO O L E U
double Det(matriz M){
    double det = 1;
    matriz** L, U;
    DecomposicaoLU(M,&L,&U);

    for(int i = 0; i < U.numColunas ;i++)
        det = det*U.elemento[i][i];
    return det;
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
            elemento[i][j] = pow((-1),(i+j))*Det(*cofij);
        }
    }
    cof->elemento = elemento;
    return cof;
}
//TA OKK
matriz* Transposta(matriz M){
    matriz* matrizT;
    double** elemento = malloc(M.numColunas* sizeof(double*));
    for(int i = 0;i < M.numColunas; i++)
        elemento[i] = malloc(M.numLinhas*sizeof(double));

    for (int i = 0; i < M.numLinhas; i++) {
        for (int j = i+1; j < M.numColunas; j++) {
            if (j != i) {
                elemento[i][j] = M.elemento[j][i];
            }
        }
    }
    matrizT->elemento = elemento;
    return matrizT;
}

matriz* Inversa(matriz M){

    return multiplicaConstante(Det(M),*Transposta(*(Cofatores(M))));
}
