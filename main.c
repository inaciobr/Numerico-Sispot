#include <stdio.h>
#include "node.h"
#include "barra.h"
#include "matriz.h"

int main() {
    matriz M = criaMatriz(3, 3);
    matriz v = criaMatriz(3,1);
    matriz kM;
    matriz LU;

    double test[3][3] = { {1.0, 2.0, 1.0},
                          {1.0, 1.0, 2.0},
                          {2.0, 1.0, 1.0} };
    double vec[3][1] = { {2.0},{3.0},{1.0} };

    for (int i = 0; i < M.numLinhas; i ++)
        for (int j = 0; j <  M.numColunas; j++)
            M.elemento[i][j] = test[i][j];
    for (int i = 0; i < v.numLinhas; i ++)
        for (int j = 0; j <  v.numColunas; j++)
            v.elemento[i][j] = vec[i][j];
    LU = multVetor(M,v);

 /**
    int *P;
    LU = decomposicaoLU(M, &P);

    kM = inversa(M);
**/
    printf ("%d, %d\n", LU.numLinhas, LU.numColunas);
    for (int i = 0; i < LU.numLinhas; i ++) {
        for (int j = 0; j <  LU.numColunas; j++)
            printf("%10.4lf ", LU.elemento[i][j]);
        printf("\n");
    }
 /**
    for (int i = 0; i < LU.numColunas; i++)
        printf("- %d\n", P[i]);

    printf("%lf", det(M));


**/
    return 0;
}
