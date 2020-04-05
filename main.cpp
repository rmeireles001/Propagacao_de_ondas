#include <iostream>
#include <cstdio>
#include "propagacao.h"
#include "aco.h"

using namespace std;

/*int main(){
	int n = fim-inicio+3;
	int laco;
	propagacao p;
	srand(time(NULL));
	p.dados_experimentais("eco.txt", n, inicio);
	p.estimativa_inicial(n);
	aco a;
	a.get_data(0,1,0.05);
	//a.run(2, &p);
	laco = 3;
	while(laco<=n){
		p.config_area(laco);
		a.run(laco, &p);
		p.attr_config(laco, a.get_var());
		p.atualizar_area(laco);
		cout << p.posicao[laco] << "   " << p.config[laco] << endl;
		laco++;
	}
	a.end();
	for(int i=3; i<=n; i++){
		cout << p.posicao[i] << "   " << p.config[i] << endl;
	}
	return 0;
}*/

/*int main(){
	clock_t start, end;
	aco a;
	propagacao p;
	srand(time(NULL));
	p.inserir("areas.txt");
	p.prob_direto(1000, p.A, p.Gexp);
	p.prob_direto(480, p.A, p.G);
	a.get_data(0,1,0.05);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		a.run(i, &p);
		p.attr_config(i, a.get_var());
		p.prob_inverso(i, p.config);
		cout << p.posicao[i] << "   " << p.config[i] << endl;
	}
	end = clock();
	FILE *resultados = fopen("areas_inverso.txt", "w");
	for(int i=inicio; i<=fim; i++){
		fprintf(resultados, "%d\t\t%.15e\n", (int) p.posicao[i], p.config[i]);
	}
	fprintf(resultados, "Tempo gasto: %lf\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
	fclose(resultados);
	a.end();
}*/

void teste(int inicio, int fim, string m, string range){
	propagacao p;
	p.run_lj("areas.txt", "areas_lj#"+m+"R"+range+".txt", inicio, fim);
	clock_t start, end;
	propagacao p2;
	p2.inserir("areas.txt");
	p2.prob_direto(1000, p2.Gexp);
	p2.prob_direto(inicio, p2.G);
	aco a;
	a.get_data(0, 1, 0.05);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		a.run(i, &p2);
		p2.atribuirA(i, a.get_var());
		p2.prob_inverso(i);
	}
	end = clock();
	p2.escrever_txt("areas_aco#"+m+"R"+range+".txt", "ANTCOLONY OPTIMIZATION", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	a.end();
}

int main(){
	srand(time(NULL));
	teste(1,100, "1", "100");
	teste(455,554, "2", "100");
	teste(991,1000, "3", "100");

	teste(1,250, "1", "250");
	teste(381,630, "2", "250");
	teste(751,1000, "3", "250");

	teste(1,500, "1", "500");
	teste(256,755, "2", "500");
	teste(501,1000, "3", "500");

	teste(1,1000, "##", "1000");
}