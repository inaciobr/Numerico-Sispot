#include <stdio.h>

#include "rede.h"
#include "numerico.h"

int main() {
	char *redes[4] = { "Redes/1_Stevenson/1_Stevenson",
					   "Redes/2_Reticulada/2_Reticulada",
		               "Redes/3_Distribuicao_Primaria/3_Distribuicao_Primaria",
		               "Redes/4_Distribuicao_Pri_Sec/4_Distribuicao_Primaria_Secundaria" };

	int num = 0;

	printf("EP 1 - Fluxo de potencia em redes eletricas pelo metodo de Newton\n");
	printf("Analise de: %s\n\n", redes[num]);

    rede *r;
    r = leituraRede(redes[num]);

    fluxoDePotenciaNewton(r);
	printDadosRede(r);

	system("pause");

    return 0;
}
