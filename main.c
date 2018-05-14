#include <stdio.h>

#include "node.h"
#include "barra.h"
#include "matriz.h"
#include "numerico.h"

void Fx(double x[3], double r[3]) {
    r[0] = 3*x[0] - cos(x[1]*x[2]) - 1/2.;
    r[1] = x[0]*x[0] - 81*(x[1] + 0.1)*(x[1] + 0.1) + sin(x[2]) + 1.06;
    r[2] = exp(-x[0]*x[1]) + 20*x[2] + (10*PI - 3)/3;
}

void dF1(double x[3], double r[3]) {
    r[0] = 3;
    r[1] = x[2]*sin(x[1]*x[2]);
    r[2] = x[1]*sin(x[1]*x[2]);
}

void dF2(double x[3], double r[3]) {
    r[0] = 2*x[0];
    r[1] = -162*(x[1] + 0.1);
    r[2] = cos(x[2]);
}

void dF3(double x[3], double r[3]) {
    r[0] = -x[1]*exp(-x[0]*x[1]);
    r[1] = -x[0]*exp(-x[0]*x[1]);
    r[2] = 20;
}

int main() {

    for (int i = 0; i < 1000; i++) {
        double x[3] = {0.1, 0.1, -0.1};
        void (*dF[3]) = {&dF1, &dF2, &dF3};

        zeroNewton(3, x, &Fx, 3, dF);
    }

    return 0;
}
