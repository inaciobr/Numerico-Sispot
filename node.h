#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int posJ;
    int posK;
    double condutancia;
    double susceptancia;
} node;

int leituraNode(node** nodes, char arquivo[]);
