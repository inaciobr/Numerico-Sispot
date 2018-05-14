#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double **condutancia;
    double **susceptancia;
} node;

void leituraNode(node *nodes, int tamanhoMatriz, char arquivo[]);


#endif // NODE_H_INCLUDED
