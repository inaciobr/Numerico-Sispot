#ifndef NODAL_H_INCLUDED
#define NODAL_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include "matriz.h"

typedef struct {
    matriz condutancia;
    matriz susceptancia;
} nodal;

void leituraNodal(nodal *mNodal, int tamanhoMatriz, char arquivo[]);
void freeNodal(nodal *mNodal);


#endif // NODAL_H_INCLUDED
