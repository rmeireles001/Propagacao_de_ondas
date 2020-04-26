#include <stdio.h>

#define tam 1000

void insertionSort(double arr[], int n){
	int i, j, key;
	for(i=1; i<n; i++){
		key = arr[i];
		j = i - 1;
		while(j>=0&&arr[j]>key){
			arr[j+1] = arr[j];
			j = j - 1;
		}
		arr[j+1] = key;
	}
}

int main(){
	double aco[tam], lj[tam];
	int i;
	double acopior=10e-10, acomelhor=10e10, ljpior=10e-10, ljmelhor=10e10;
	int acomp=0, acopp=0, ljmp=0, ljpp=0;
	FILE *f = fopen("tempos.txt", "r");
	for(i=0; i<tam; i++){
		fscanf(f, "%lf", &aco[i]);
	}
	for(i=0; i<tam; i++){
		fscanf(f, "%lf", &lj[i]);
	}
	fclose(f);
	for(i=0; i<tam; i++){
		if(aco[i]>acopior){
			acopior = aco[i];
			acopp = i+1;
		}
		if(aco[i]<acomelhor){
			acomelhor = aco[i];
			acomp = i+1;
		}
		if(lj[i]>ljpior){
			ljpior = lj[i];
			ljpp = i+1;
		}
		if(lj[i]<ljmelhor){
			ljmelhor = lj[i];
			ljmp = i+1;
		}
	}
	insertionSort(aco, tam);
	insertionSort(lj, tam);
	FILE *s = fopen("ord.txt", "w");
	fprintf(s, "ANTCOLONY OPTIMIZATION:\nPior execução: #%d\nTempo da pior execução: %lf\nMelhor execução: #%d\nTempo da melhor execução: %lf\n", acopp, acopior, acomp, acomelhor);
	fprintf(s, "LUUS JAAKOLA:\nPior execução: #%d\nTempo da pior execução: %lf\nMelhor execução: #%d\nTempo da melhor execução: %lf\n", ljpp, ljpior, ljmp, ljmelhor);
	fprintf(s, "\n\nANTCOLONY OPTIMIZATION\n\n");
	for(i=0; i<tam; i++){
		fprintf(s, "%lf\n", aco[i]);
	}
	fprintf(s, "\nLUUS JAAKOLA\n\n");
	for(i=0; i<tam; i++){
		fprintf(s, "%lf\n", lj[i]);
	}
	fclose(s);
	return 0;
}
