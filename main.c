#include <stdio.h>
#include "node.h"
#include "barra.h"
#include "matriz.h"
#include "numerico.h"

double* dF1(double* x) {
    double *r = malloc(3*sizeof(double));

    r[0] = 3;
    r[1] = x[2]*sin(x[1]*x[2]);
    r[2] = x[1]*sin(x[1]*x[2]);

    return r;
}

double* dF2(double* x) {
    double *r = malloc(3*sizeof(double));

    r[0] = 2*x[0];
    r[1] = -162*(x[1] + 0.1);
    r[2] = cos(x[2]);

    return r;
}

double* dF3(double* x) {
    double *r = malloc(3*sizeof(double));

    r[0] = -x[1]*exp(-x[0]*x[1]);
    r[1] = -x[0]*exp(-x[0]*x[1]);
    r[2] = 20;

    return r;
}

int main() {
    double x[3] = {0.1, 0.1, -0.1};

    matriz M = matrizJacobiana(x, &dF1, &dF2, &dF3);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            printf("%15.8lf ", M.elemento[i][j]);
        printf("\n");
    }

    return 0;
}
