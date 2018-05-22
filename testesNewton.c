/** TESTES INICIAIS */

/** Teste 1*/
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

/** Teste 2*/
void F2x(double x[4], double r[4]) {
    r[0] = 4*x[0] - x[1] + x[2] - x[0]*x[3];
    r[1] = -x[0] + 3*x[1] - 2*x[2] - x[1]*x[3];
    r[2] = x[0] - 2*x[1] + 3*x[2] - x[2]*x[3];
    r[3] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1;
}

void dF2dx1(double x[4], double r[4]) {
    r[0] = 4 - x[3];
    r[1] = -1;
    r[2] = 1;
    r[3] = 2*x[0];
}

void dF2dx2(double x[4], double r[4]) {
    r[0] = -1;
    r[1] = 3 - x[3];
    r[2] = -2;
    r[3] = 2*x[1];
}

void dF2dx3(double x[4], double r[4]) {
    r[0] = 1;
    r[1] = -2;
    r[2] = 3 - x[3];
    r[3] = 2*x[2];
}

void dF2dx4(double x[4], double r[4]) {
    r[0] = -x[0];
    r[1] = -x[1];
    r[2] = -x[2];
    r[3] = 0;
}

/** Teste 3*/
void F3x20(double x[20], double r[20]) {
    for (int i = 0; i < 20; i++) {
        r[i] = i >= 1 ? -x[i - 1] : 0;
        r[i] += 2*x[i] - x[i + 1] - exp(x[i])/(20*20);
    }
}

void dF3dx(int i, double x[20], double r[20]) {
    r[i - 1] = -1;
    r[i] = 2 - exp(x[i])/(20*20);
    r[i + 1] = -1;
}

void* dF3dxx(double x[20], double r[20]) {
    void (*dF2[20]);
    for (int i = 0; i < 20; i++)
        dF2[i] = &dF3dx;

    return dF2;
}

void testesZeroNewton() {
    int iteracoes;
    printf("Testes iniciais:\n\n");

    /** Teste 1*/
    printf("1 - Minimo da funcao F(x, y) = (x - 2)^2 + (y - 3)^2\n");

    double x1[2] = {0.0, 0.0};
    void (*dF1[2]) = {&dF1dx1, &dF1dx2};

    iteracoes = zeroNewton(2, x1, &F1x, 2, &dF1);

    printf("x = %.10f, y = %.10f\n", x1[0], x1[1]);
    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);


    /** Teste 2*/
    printf("1 - Raiz da funcao o F(x1, x2, x3, x4) = (4x1 - x2  + x3 - x1*x4, -x1 + 3x2 - 2x3 - x2*x4, x1 - 2x2 + 3x3 - x3*x4, x1^2 + x2^2 + x3^2 - 1)\n");

    double x2[4] = {1.0, 1.0, 1.0, 1.0};
    void (*dF2[4]) = {&dF2dx1, &dF2dx2, &dF2dx3, &dF2dx4};

    iteracoes = zeroNewton(4, x2, &F2x, 4, &dF2);

    printf("x1 = %.10f, x2 = %.10f, x3 = %.10f, x4 = %.10f\n", x2[0], x2[1], x2[2], x2[3]);
    printf("Resultado em %d iteracoes.\n\n", iteracoes + 1);


    /** Teste 3*/
    double x3[20] = {};
}
