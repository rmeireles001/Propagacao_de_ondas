#include <iostream>
#include <cstdio>
#define runs 1000

using namespace std;

class block{
    public:
    double times[runs];
    int run[runs];
    string algorithm;
    void read(string name);
    void bubble_sort(int n);
    void swap(double *a, double *b);
    void swap(int *a, int *b);
    void final_results(FILE *f);
};

int main(){
    block aco, lj, cgrasp;
    aco.algorithm = "ANT COLONY OPTIMIZATION";
    lj.algorithm = "LUUS-JAAKOLA";
    cgrasp.algorithm = "C-GRASP";
    FILE *output = fopen("output.txt", "w");
    fprintf(output, "Beam#1\n-------------------\n\n");
    aco.read("aco#1.txt");
    lj.read("lj#1.txt");
    cgrasp.read("cgrasp#1.txt");
    aco.bubble_sort(runs);
    aco.final_results(output);
    lj.bubble_sort(runs);
    lj.final_results(output);
    cgrasp.bubble_sort(runs);
    cgrasp.final_results(output);
    aco.read("aco#2.txt");
    lj.read("lj#2.txt");
    cgrasp.read("cgrasp#2.txt");
    fprintf(output, "Beam#2\n-------------------\n\n");
    aco.bubble_sort(runs);
    aco.final_results(output);
    lj.bubble_sort(runs);
    lj.final_results(output);
    cgrasp.bubble_sort(runs);
    cgrasp.final_results(output);
    fclose(output);
    return 0;
}

void block::read(string name){
    FILE *f = fopen(name.c_str(), "r");
    for(int i=0; i<runs; i++){
        fscanf(f, "%lf", &times[i]);
        run[i] = i;
    }
    fclose(f);
}

void block::bubble_sort(int n){
    if(n<1)
        return;
    for(int i=0; i<n-1; i++){
        if(times[i]>times[i+1]){
            swap(&times[i], &times[i+1]);
            swap(&run[i], &run[i+1]);
        }
    }
    bubble_sort(n-1);
}
void block::swap(double *a, double *b){
    *a = *a+*b;
    *b = *a-*b;
    *a = *a-*b;
}

void block::swap(int *a, int *b){
    *a = *a+*b;
    *b = *a-*b;
    *a = *a-*b;
}

void block::final_results(FILE *f){
    fprintf(f, "Algorithm: %s\n", algorithm.c_str());
    fprintf(f, "Worst run: %d\n", run[runs-1]);
    fprintf(f, "Worst time: %lf\n", times[runs-1]);
    fprintf(f, "Best run: %d\n", run[0]);
    fprintf(f, "Best time: %lf\n\n", times[0]);
}