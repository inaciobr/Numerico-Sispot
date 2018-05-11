#include <stdio.h>

#include "node.h"
#include "barra.h"
#include "matriz.h"
#include "numerico.h"

double *Fx(double x[3]) {
    double *r = malloc(3 * sizeof(double));

    r[0] = 3*x[0] - cos(x[1]*x[2]) - 1/2.;
    r[1] = x[0]*x[0] - 81*(x[1] + 0.1)*(x[1] + 0.1) + sin(x[2]) + 1.06;
    r[2] = exp(-x[0]*x[1]) + 20*x[2] + (10*M_PI - 3)/3;

    return r;
}

double* dF1(double x[3]) {
    double *r = malloc(3 * sizeof(double));

    r[0] = 3;
    r[1] = x[2]*sin(x[1]*x[2]);
    r[2] = x[1]*sin(x[1]*x[2]);

    return r;
}

double* dF2(double x[3]) {
    double *r = malloc(3 * sizeof(double));

    r[0] = 2*x[0];
    r[1] = -162*(x[1] + 0.1);
    r[2] = cos(x[2]);

    return r;
}

double* dF3(double x[3]) {
    double *r = malloc(3 * sizeof(double));

    r[0] = -x[1]*exp(-x[0]*x[1]);
    r[1] = -x[0]*exp(-x[0]*x[1]);
    r[2] = 20;

    return r;
}

int main() {
    double x[3] = {0.1, 0.1, -0.1};
    matriz dF = matrizFuncao(3, x, 3, &dF1, &dF2, &dF3);
    matriz F = matrizFuncao(3, x, 1, &Fx);

    for (int i = 0; i < dF.numLinhas; i++) {
        for (int j = 0; j < dF.numColunas; j++)
            printf("%15.8lf ", dF.elemento[i][j]);

        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < F.numLinhas; i++) {
        for (int j = 0; j < F.numColunas; j++)
            printf("%15.8lf ", F.elemento[i][j]);

        printf("\n");
    }

    return 0;
}
