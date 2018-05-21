#include <stdio.h>

#include "rede.h"
#include "numerico.h"

int main() {
    int numRedes = 4;
	char *redes[4] = { "Redes/1_Stevenson/1_Stevenson",
					   "Redes/2_Reticulada/2_Reticulada",
		               "Redes/3_Distribuicao_Primaria/3_Distribuicao_Primaria",
		               "Redes/4_Distribuicao_Pri_Sec/4_Distribuicao_Primaria_Secundaria" };

	printf("EP 1 - Fluxo de potencia em redes eletricas pelo metodo de Newton\n");

	while (1) {
        printf("Digite o numero referente ao que deseja fazer:\n");

        printf("0 - Sair\n");
        for (int i = 0; i < numRedes; i++) {
            printf("%d - %s\n", i + 1, redes[i]);
        }

        int menu;
        scanf("%d", &menu);

        if (menu == 0)
            return 0;

        rede *r;
        r = leituraRede(redes[menu - 1]);

        fluxoDePotenciaNewton(r);
        printDadosRede(r);
        arquivarDadosRede(r);

        freeRede(r);
	}


    return 0;
}
