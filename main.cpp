#include <iostream>
#include <cstdio>
#include "propagacao.h"
#include "aco.h"

using namespace std;

void teste(int inicio, int fim, string m, string range){
	//propagacao p;
	//p.run_lj("areas.txt", "areas_lj#"+m+"R"+range+".txt", inicio, fim);
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
	/*teste(1,100, "1", "100");
	teste(455,554, "2", "100");
	teste(991,1000, "3", "100");

	teste(1,250, "1", "250");
	teste(381,630, "2", "250");
	teste(751,1000, "3", "250");

	teste(1,500, "1", "500");
	teste(256,755, "2", "500");
	teste(501,1000, "3", "500");*/

	teste(1,1000, "##", "1000");
}
