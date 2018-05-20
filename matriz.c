#include "matriz.h"

/**
 * Aloca o espaço necessário para uma matriz de acordo com o
 * número de linhas e de colunas  e de colunas solicitado.
 * Todos os elementos da matriz são do tipo double.
 */
matriz criaMatriz(int numLinhas, int numColunas) {
    matriz M;
    M.numLinhas = numLinhas;
    M.numColunas = numColunas;

    M.elemento = malloc(numLinhas * sizeof(double *));

    for (int i = 0; i < numLinhas; i++)
        M.elemento[i] = calloc(numColunas, sizeof(double));

    return M;
}

/**
 * Desaloca a matriz alocada dinamicamente.
 */
void freeMatriz(matriz *A) {
    for (int i = 0; i < A->numLinhas; i++)
        free(A->elemento[i]);

    free(A->elemento);

    A->numLinhas = 0;
    A->numColunas = 0;
}

/**
 * Multiplica todos os elementos de uma matriz por uma constante.
 * A matriz utilizada é alterada no interior da função.
 */
matriz multiplicaConstante(matriz A, double constante) {
    for (int i = 0; i < A.numLinhas; i++)
        for (int j = 0; j < A.numColunas; j++)
            A.elemento[i][j] *= constante;

    return A;
}

/**
 * Calcula o determinante de uma matriz e o retorna com o tipo double.
 * O determinante é calculado a partir do produto da diagonal
 * principal da matriz após a decomposição LU.
 */
double det(matriz M) {
    if (M.numColunas != M.numLinhas)
        return 0;

    int *P;
    matriz LU = copiaMatriz(M);
    LU = decomposicaoLU(LU, &P);

    int N = P[M.numLinhas];
    free(P);

    double det = 1;
    for (int i = 0; i < M.numLinhas; i++)
        det *= LU.elemento[i][i];

    freeMatriz(&LU);

    return det * (N%2 ? -1 : 1);
}

/**
 * Calcula o determinante de uma matriz e o retorna com o tipo double.
 * O determinante é calculado a partir do produto da diagonal
 * principal da matriz após a decomposição LU.
 */
matriz transposta(matriz M) {
    matriz matrizT = criaMatriz(M.numColunas, M.numLinhas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
                matrizT.elemento[j][i] = M.elemento[i][j];

    freeMatriz(&M);

    return matrizT;
}

/**
 * Cria uma nova matriz igual à passada por argumento, porém
 * realizando a aloação novamente, para casos em que a matriz
 * original é alterada, mas se deseja manter uma cópia.
 */
matriz copiaMatriz(matriz M) {
    matriz copiaM = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            copiaM.elemento[i][j] = M.elemento[i][j];

    return copiaM;
}

/**
 * Realiza a decomposição LU com pivotação de uma matriz.
 * O vetor de inteiros permutacoes passado como argumento
 * armazena o resultado dessas permutações e sua última linha
 * informa o número de permutações realizadas.
 * Há apenas uma matriz resultante onde os elementos abaixo
 * da diagonal principal são referentes à matriz L e os
 * elementos da diagonal principal e acima dela são referentes
 * à matriz U.
 */
matriz decomposicaoLU(matriz M, int **permutacoes) {
    matriz LU = M;
    int *P = calloc((M.numLinhas + 1), sizeof(double));
	double maxPivo, abs;
	double *temp;

    for (int k = 0; k < M.numLinhas; k++) {
        maxPivo = -1.0;
        abs = 0.0;

        for (int i = k; i < M.numLinhas; i++) {
            for (int j = 0; j < k; j++)
                LU.elemento[i][k] -= LU.elemento[i][j] * LU.elemento[j][k];

            /* Salva o valor em módulo do valor atual. */
            abs = fabs(LU.elemento[i][k]);

            /* Salva o valor do maior pivô até o momento */
            if (abs > maxPivo) {
                maxPivo = abs;
                P[k] = i;
            }
        }

        /* Troca duas linhas da matriz. */
        if (P[k] != k) {
            temp = LU.elemento[k];
            LU.elemento[k] = LU.elemento[P[k]];
            LU.elemento[P[k]] = temp;

            P[M.numLinhas]++;
        }

        for (int j = k + 1; j < M.numLinhas; j++) {
            for (int i = 0; i < k; i++)
                LU.elemento[k][j] -= LU.elemento[k][i] * LU.elemento[i][j];

            LU.elemento[j][k] /= LU.elemento[k][k];
        }
    }

    *permutacoes = P;
    return LU;
}

/**
 * Realiza a permutação de elementos de uma matriz
 * seguindo o vetor criado pela função de decomposição LU.
 * O vetor é lido da primeira linha para a última, caso o
 * valor informado em uma linha seja diferente do número
 * desta linha, as linhas referidas são trocadas.
 */
void permutaLinhasMatriz(matriz A, int P[]) {
    for (int i = 0; i < A.numLinhas; i++)
        if (i != P[i]) {
            double *temp = A.elemento[i];
            A.elemento[i] = A.elemento[P[i]];
            A.elemento[P[i]] = temp;
        }
}

/**
 * Calcula o valor do cofator de um determinado elemento
 * de uma matriz.
 */
double cofator(matriz M, int linha, int coluna) {
    matriz mCofatores = criaMatriz(M.numLinhas - 1, M.numColunas - 1);
    double resultado;

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            if (i != linha && j != coluna)
                mCofatores.elemento[(i <= linha ? i : i - 1)][(j <= coluna ? j : j - 1)] = M.elemento[i][j];

    resultado = ((linha + coluna)%2 ? -1 : 1) * det(mCofatores);

    freeMatriz(&mCofatores);

    return resultado;
}

/**
 * Cria uma matriz com todos os cofatores de uma
 * determinada matriz.
 * A matriz original não é alterada.
 */
matriz matrizCofatores(matriz M) {
    matriz cofatores = criaMatriz(M.numLinhas, M.numColunas);

    for (int i = 0; i < M.numLinhas; i++)
        for (int j = 0; j < M.numColunas; j++)
            cofatores.elemento[i][j] = cofator(M, i, j);

    return cofatores;
}

/**
 * Calcula a inversa de uma matriz.
 * A matriz original é não alterada.
 */
matriz inversa(matriz M) {
    if (M.numColunas != M.numLinhas)
        return M;

    matriz cof = matrizCofatores(M);
    matriz inversa = multiplicaConstante(transposta(cof), 1/det(M));

    return inversa;
}

/**
 * Realiza o produto matricial entre duas matrizes.
 * As matrizes não são alteradas.
 */
matriz produtoMatriz(matriz M1, matriz M2){
    matriz produto = criaMatriz(M1.numLinhas, M2.numColunas);

    for (int i = 0; i < M1.numLinhas; i++)
        for (int j = 0; j < M2.numColunas; j++)
                for (int k = 0; k < M1.numColunas; k++)
                    produto.elemento[i][j] += M1.elemento[i][k] * M2.elemento[k][j];

    return produto;
}

/**
 * Resolve sistema da forma Ax = b utilizando a
 * decomposição LU.
 * A matriz A é alterada.
 */
matriz resolveSistemaLinear(matriz A, matriz b) {
    int *P;
    double soma;
    matriz LU = decomposicaoLU(A, &P);

    /* Permuta o vetor b para seguir a ordem da matriz A decomposta */
    permutaLinhasMatriz(b, P);

    free(P);

    /* Resolve primeiramente o sistema Ly = b */
    matriz y = criaMatriz(A.numLinhas, 1);
    for (int i = 0; i < A.numLinhas; i++) {
        soma = 0.0;

        for (int j = 0; j < i; j++)
            soma += LU.elemento[i][j] * y.elemento[j][0];

        y.elemento[i][0] = (b.elemento[i][0] - soma);
    }

    /* Resolve o sistema Ux = y */
    matriz x = criaMatriz(A.numLinhas, 1);
    for (int i = A.numLinhas - 1; i >= 0; i--) {
        soma = 0.0;

        for (int j = i; j < A.numColunas; j++)
            soma += LU.elemento[i][j] * x.elemento[j][0];

        x.elemento[i][0] = (y.elemento[i][0] - soma) / LU.elemento[i][i];
    }

    freeMatriz(&y);

    return x;
}

/**
 * Exibe todos os elementos de uma matriz.
 */
void printMatriz(matriz M) {
    for (int i = 0; i < M.numLinhas; i++) {
        for (int j = 0; j < M.numColunas; j++)
            printf("%15.2f ", M.elemento[i][j]);

        printf("\n");
    }
}

