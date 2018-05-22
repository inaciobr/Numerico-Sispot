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
        printf("%d - Testar metodo de Newton\n", numRedes + 1);

        int menu;
        scanf("%d", &menu);

        if(menu == numRedes + 1) {
            testesZeroNewton();
            continue;
        }
        else if (menu == 0 || menu > numRedes + 1)
            return 0;

        rede *r;
        r = leituraRede(redes[menu - 1]);

        int iteracoes;
        iteracoes = fluxoDePotenciaNewton(r);
        printf("\nResultado em %d iteracoes.\n", iteracoes + 1);

        printDadosRede(r);
        arquivarDadosRede(r);

        freeRede(r);
	}


    return 0;
}
