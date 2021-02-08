#include <iostream>
#include <cstdio>
#include "propagacao.h"
#include "aco.h"

#define contagem 10

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

double run_aco(int inicio, int fim, string arq_areas, string runn){
	clock_t start, end;
	propagacao p;
	p.inserir(arq_areas);
	p.prob_direto(1000, p.Gexp);
	p.prob_direto(inicio, p.G);
	aco a;
	a.get_data(0, 1, 0.05);
	start = clock();
	cout << "Foi\n\n";
	for(int i=inicio; i<=fim; i++){
		a.run(i, &p);
		p.atribuirA(i, a.get_var());
		p.prob_inverso(i);
	}
	end = clock();
	p.escrever_txt(runn, "ANTCOLONY OPTIMIZATION", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	a.end();
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

double run_lj(int inicio, int fim, string arq_areas, string runn){
	clock_t start, end;
	double img;
	propagacao p;
	p.inserir(arq_areas);
	p.prob_direto(1000, p.Gexp);
	p.prob_direto(inicio, p.G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		img = p.luus_jaakola(i);
		p.atribuirA(i, img);
		p.prob_inverso(i);
	}
	end = clock();
	p.escrever_txt(runn, "LUUS JAAKOLA", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

double run_cgrasp(int inicio, int fim, string arq_areas, string runn){
	clock_t start, end;
	double img;
	propagacao p;
	p.inserir(arq_areas);
	p.prob_direto(1000, p.Gexp);
	p.prob_direto(inicio, p.G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		img = p.cgrasp(0, 1, 0.05, i);
		p.atribuirA(i, img);
		p.prob_inverso(i);
	}
	end = clock();
	p.escrever_txt(runn, "C-GRASP", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}


void count_aco(FILE *resposta, string arq_areas, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "ANTCOLONY OPTIMIZATION ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=1; i<=contagem; i++){
		cout << "entou\n";
		tempos[i] = run_aco(inicio, fim, arq_areas, "aco_run#"+int2str(i)+".txt");
		fprintf(resposta, "%lf\t", tempos[i]);
		cout << tempos[i] << endl;
	}
	cout << "fora do for" << endl;
	for(int i=0; i<contagem; i++){
		med = med + tempos[i];
	}
	med=med/contagem;
	for(int i=0; i<contagem; i++){
		var = var + pow(tempos[i] - med, 2);
	}
	var = var/(contagem);
	cout << "Está para escrever" << endl;
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
	cout << "Fim dos testes ACO\n";
}

void count_lj(FILE *resposta, string arq_areas, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "LUUS JAAKOLA ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=1; i<=contagem; i++){
		tempos[i] = run_lj(inicio, fim, arq_areas, "lj_run#"+int2str(i)+".txt");
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
	var = var/(contagem);
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
}

void count_cgrasp(FILE *resposta, string arq_areas, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "C-GRASP ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=1; i<=contagem; i++){
		tempos[i] = run_cgrasp(inicio, fim, arq_areas, "cgrasp_run#"+int2str(i)+".txt");
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
	var = var/(contagem);
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

double run_ord(int inicio, int fim, string arq_areas, string runn){
	clock_t start, end;
	double img;
	propagacao p;
	p.inserir(arq_areas);
	p.prob_direto(1000, p.Gexp);
	p.prob_direto(inicio, p.G);
	p.prob_inverso();
	ord *vetor = new ord[21];
	start = clock();
	for(int i=inicio; i<=fim; i++){
		cout << i << endl;
		img = p.busca_ordenada(vetor, i, 0.05, 0, 1);
		p.atribuirA(i, img);
		p.prob_inverso(i);
		if(vetor[0].ggexp>1e-3){
			cout << "ls" << endl;
			double llimit = img-0.05;
			double ulimit = img+0.05;
			if(llimit<0.00){
				llimit = 0.00;
			}
			if(ulimit>1.00){
				ulimit = 1.00;
			}
			img = p.busca_ordenada(vetor, i, 0.01, llimit, ulimit);
			p.atribuirA(i, img);
			p.prob_inverso(i);
		}
	}
	end = clock();
	p.escrever_txt(runn, "ORD_SEARCH", (double)(end-start)/(double)(CLOCKS_PER_SEC), 1, fim);
	//p.imprime_area(arv, 1, 0.05, 498);
	delete vetor;
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

void count_ord(FILE *resposta, string arq_areas, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "ORD ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=1; i<=contagem; i++){
		tempos[i] = run_ord(inicio, fim, arq_areas, "ords_run#"+int2str(i)+".txt");
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
	var = var/(contagem);
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
}

int main(){
	srand(time(NULL));
	FILE *results = fopen("res.txt", "w");
	count_ord(results, "areas.txt", 1, 1000);
	fclose(results);
	system("mkdir results\nmv ords_run#*.txt results\nmv res.txt results");
	results = fopen("res.txt", "w");
	count_ord(results, "areas2.txt", 1, 1000);
	fclose(results);
	system("mkdir results2\nmv ords_run#*.txt results2\nmv res.txt results2");
	return 0;
}