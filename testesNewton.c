#include "testesNewton.h"

/** TESTES INICIAIS */

/** Teste 1*/
void F1(matriz *M, double x[2]) {
    M->elemento[0][0] = 2*(x[0] - 2);
    M->elemento[1][0]= 2*(x[1] - 3);
}

void JF1(matriz *M, double x[2]) {
    /* dF1 / dx */
    M->elemento[0][0] = 2;
    M->elemento[1][0] = 0;

    /* dF1 / dy */
    M->elemento[0][1] = 0;
    M->elemento[1][1] = 2;
}

/** Teste 2*/
void F2(matriz *M, double x[4]) {
    M->elemento[0][0] = 4*x[0] - x[1] + x[2] - x[0]*x[3];
    M->elemento[1][0] = -x[0] + 3*x[1] - 2*x[2] - x[1]*x[3];
    M->elemento[2][0] = x[0] - 2*x[1] + 3*x[2] - x[2]*x[3];
    M->elemento[3][0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1;
}

void JF2(matriz *M, double x[4]) {
    /* dF2 / dx1 */
    M->elemento[0][0] = 4 - x[3];
    M->elemento[1][0] = -1;
    M->elemento[2][0] = 1;
    M->elemento[3][0] = 2*x[0];

    /* dF2 / dx2 */
    M->elemento[0][1] = -1;
    M->elemento[1][1] = 3 - x[3];
    M->elemento[2][1] = -2;
    M->elemento[3][1] = 2*x[1];;

    /* dF2 / dx3 */
    M->elemento[0][2] = 1;
    M->elemento[1][2] = -2;
    M->elemento[2][2] = 3 - x[3];
    M->elemento[3][2] = 2*x[2];

    /* dF2 / dx4 */
    M->elemento[0][3] = -x[0];
    M->elemento[1][3] = -x[1];
    M->elemento[2][3] = -x[2];
    M->elemento[3][3] = 0;
}

/** Teste 3*/
void F3(matriz *M, double x[]) {
    int N = M->numLinhas;

    for (int i = 0; i < N - 1; i++)
        M->elemento[0][0] = (i > 1 ? -x[i - 1] : 0.0) + 2.0*x[i] - x[i + 1] - exp(x[i])/(N*N);

}

void JF3(matriz *M, double x[]) {
    int N = M->numLinhas;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (i == j)
                M->elemento[i][j] = 2.0 - exp(x[i])/(N*N);
            else if (i == j + 1 || i == j - 1)
                M->elemento[i][j] = -1.0;
            else
                M->elemento[i][j] = 0.0;
}
void testesZeroNewton() {
    int iteracoes;
    printf("Testes iniciais:\n\n");

    /** Teste 1*/
    printf("1 - Minimo da funcao F(x, y) = (x - 2)^2 + (y - 3)^2\n");

    double x1[2] = {1.0, 1.0};
    iteracoes = zeroNewton(2, x1, &F1, &JF1);

    printf("x = %.10f, y = %.10f\n", x1[0], x1[1]);
    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);


    /** Teste 2*/
    printf("2 - Raiz da funcao o F(x1, x2, x3, x4) = (4x1 - x2  + x3 - x1*x4, -x1 + 3x2 - 2x3 - x2*x4, x1 - 2x2 + 3x3 - x3*x4, x1^2 + x2^2 + x3^2 - 1)\n");

    double x2[4] = {1.0, 1.0, 1.0, 1.0};
    iteracoes = zeroNewton(4, x2, &F2, &JF2);

    printf("x1 = %.10f, x2 = %.10f, x3 = %.10f, x4 = %.10f\n", x2[0], x2[1], x2[2], x2[3]);
    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);


    /** Teste 3*/
    printf("3 - Sistema n - 1 x n - 1\n");

    /** Caso N = 20 */
    printf("a) N = 20\n");
    double x3a[20] = {};
    iteracoes = zeroNewton(20, x3a, &F3, &JF3);

    printf("x = {");
    for (int i = 0; i < 20 - 1; i++)
        printf("%.3f, ", x3a[i]);
    printf("%.3f}\n", x3a[20 - 1]);

    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);

    /** Caso N = 40 */
    printf("a) N = 40\n");
    double x3b[40] = {};
    iteracoes = zeroNewton(40, x3b, &F3, &JF3);

    printf("x = {");
    for (int i = 0; i < 40 - 1; i++)
        printf("%.3f, ", x3b[i]);
    printf("%.3f}\n", x3b[40 - 1]);

    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);

    /** Caso N = 80 */
    printf("a) N = 80\n");
    double x3c[80] = {};
    iteracoes = zeroNewton(80, x3c, &F3, &JF3);

    printf("x = {");
    for (int i = 0; i < 80 - 1; i++)
        printf("%.3f, ", x3c[i]);
    printf("%.3f}\n", x3c[80 - 1]);

    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);
}
