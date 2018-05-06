#include <stdio.h>
#include "node.h"
#include "barra.h"
#include "matriz.h"

int main() {
    matriz M = criaMatriz(3, 3);
    matriz kM;
    matriz L, U;
    double test[3][3] = { {350.0, 17.0, 45.0},
                          {0.0, 18.0, 98.0},
                          {15.0, 32.0, 103.0} };

    for (int i = 0; i < M.numLinhas; i ++)
        for (int j = 0; j <  M.numColunas; j++)
            M.elemento[i][j] = test[i][j];

    decomposicaoLU(M, &L, &U);

    kM = inversa(M);

    printf ("%d, %d\n", M.numLinhas, M.numColunas);
    for (int i = 0; i < kM.numLinhas; i ++) {
        for (int j = 0; j <  kM.numColunas; j++)
            printf("%10.4lf ", kM.elemento[i][j]);
        printf("\n");
    }

    printf("%lf", det(M));



    return 0;
}
