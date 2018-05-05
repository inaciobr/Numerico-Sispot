#include "matriz.h"

//void criaMatriz(){}

void DecomposicaoLU(matriz** A,matriz** L,matriz** U){

}

double Det(matriz** matriz){
    double det = 1;
    matriz* L,U = NULL;
    DecomposicaoLU(&matriz,&L,&U);

    for(int i = 0; i < U.numColunas ;i++)
        det = det*U.elemento[i][i];
    return det;
}


matriz** Cofatores(matriz** matriz){

}

matriz** Inversa(matriz** matriz){

}
