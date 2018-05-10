#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int posJ;
    int posK;
    double condutancia;
    double susceptancia;
} node;

int leituraNode(node** nodes, char arquivo[]);


#endif // NODE_H_INCLUDED
