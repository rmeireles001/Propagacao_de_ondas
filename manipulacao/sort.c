#include <stdio.h>

#define tam 100

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
	FILE *f = fopen("tempos.txt", "r");
	for(i=0; i<tam; i++){
		fscanf(f, "%lf", aco[i]);
	}
	for(i=0; i<tam; i++){
		fscanf(f, "%lf", lj[i]);
	}
	fclose(f);
	insertionSort(aco, tam);
	insertionSort(lj, tam);
	FILE *s = fopen("ord.txt", "w");
	fprintf(s, "ANTCOLONY OPTIMIZATION\n\n");
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
