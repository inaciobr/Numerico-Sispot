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

///

void F2x(double x[4], double r[4]) {
    r[0] = 4*x[0] - x[1] + x[2] +x[0]*x[3];
    r[1] = -x[0] + 3*x[1] - 2*x[2] -x[1]*x[3];
    r[2] = x[0] - 2*x[1] + 3*x[2] - x[2]*x[3];
    r[3] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1;
}

void dF2dx1(double x[4], double r[4]) {
    r[0] = 4 + x[3];
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
    r[0] = x[0];
    r[1] = -x[1];
    r[2] = -x[2];
    r[3] = 0;
}

void testesZeroNewton() {
    double x1[4] = {1,1,1,1};
    void (*dF[4]) = { &dF2dx1, &dF2dx2, &dF2dx3, &dF2dx4};

    zeroNewton(4, x1, &F2x, 4, &dF);

    for (int i = 0; i < 3; i++)
        printf("%15.10f ", x1[i]);

    printf("\n\n");

}
