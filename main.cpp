#include <iostream>
#include <cstdio>
#include "propagacao.h"
#include "aco.h"

#define contagem 10
#define ini 450
#define fi 550

using namespace std;
string int2str(int num);

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

double run_aco(int inicio, int fim, string runn){
	clock_t start, end;
	propagacao p;
	p.inserir("areas.txt");
	p.prob_direto(1000, p.Gexp);
	p.prob_direto(inicio, p.G);
	aco a;
	a.get_data(0, 1, 0.05);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		a.run(i, &p);
		p.atribuirA(i, a.get_var());
		p.prob_inverso(i);
	}
	end = clock();
	p.escrever_txt(runn, "ANTCOLONY OPTIMIZATION", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

double run_lj(int inicio, int fim, string runn){
	clock_t start, end;
	double img;
	propagacao p;
	p.inserir("areas.txt");
	p.prob_direto(1000, p.Gexp);
	p.prob_direto(inicio, p.G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		cout << "Execução#" << runn << " area#" << i << endl;
		img = p.luus_jaakola(i);
		p.atribuirA(i, img);
		p.prob_inverso(i);
	}
	end = clock();
	p.escrever_txt(runn, "LUUS JAAKOLA", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

void count_aco(FILE *resposta, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "ANTCOLONY OPTIMIZATION ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=0; i<contagem; i++){
		tempos[i] = run_aco(inicio, fim,"aco_run#"+int2str(i)+".txt");
		fprintf(resposta, "%lf\t", tempos[i]);
		cout << tempos[i] << endl;
	}
	for(int i=0; i<contagem; i++){
		med = med + tempos[i];
	}
	med=med/contagem;
	for(int i=0; i<contagem; i++){
		var = var + pow(tempos[i] - med, 2);
	}
	var = var/(contagem-1);
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
}

void count_lj(FILE *resposta, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "LUUS JAAKOLA ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=0; i<contagem; i++){
		tempos[i] = run_lj(inicio, fim, "lj_run#"+int2str(i)+".txt");
		fprintf(resposta, "%lf\t", tempos[i]);
		cout << tempos[i] << endl;
	}
	for(int i=0; i<contagem; i++){
		med = med + tempos[i];
	}
	med=med/contagem;
	for(int i=0; i<contagem; i++){
		var = var + pow(tempos[i] - med, 2);
	}
	var = var/(contagem-1);
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
}

string int2str(int num){
	string ret="";
	int bit;
	while(num>0){
		bit = num%10;
		num = (num-bit)/10;
		ret = (char)(bit+48)+ret;
	}
	return ret;
}

int main(){
	srand(time(NULL));

	FILE *saida = fopen("saida.txt", "w");

	//count_aco(saida, 1, 1000);
	count_lj(saida, 1, 1000);
	
	fclose(saida);

}
