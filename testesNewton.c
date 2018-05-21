#include "testesNewton.h"

// TESTES INICIAIS
void F1x(double x[2], double r[2]) {
    r[0] = 2*(x[0] - 2);
    r[1] = 2*(x[1] - 3);
}

void dF1dx1(double x[2], double r[2]) {
    r[0] = 2;
    r[1] = 0;
}

void dF1dx2(double x[2], double r[2]) {
    r[0] = 0;
    r[1] = 2;
}

void testesZeroNewton() {
    double x1[2] = {0.1, -0.1};
    void (*dF1[2]) = {&dF1dx1, &dF1dx2};

    zeroNewton(2, x1, &F1x, 2, &dF1);

    for (int i = 0; i < 2; i++)
        printf("%15.10f ", x1[i]);

    printf("\n\n");

}
