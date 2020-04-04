#include "variable.h"

void variable::init(double llimit, double ulimit, double increment){
	this->ulimit = ulimit;
	this->llimit = llimit;
	this->increment = increment;
	find_nvalue();
	trail = allocate(nvalues);
	prob = allocate(nvalues);
}

void variable::free(){
	delete trail;
	delete prob;
}

void variable::find_nvalue(){
	double nv, quotient;
	int divisor;
	nv = (ulimit-llimit)/increment;
	divisor = (int) nv;
	nvalues = 1+divisor;
	quotient = nv/((double)divisor);
	if(quotient != 1.0f){
		if((fabs((quotient - 1.000000000000000) ) > 0.01)){
			nvalues++;
		}
	}
}

double *variable::allocate(int size){
	double *vet = new double[size];
	if(vet == NULL){
		cout << "Memória insuficiente" << endl;
	}
	reset(vet, size);
	return vet;
}

void variable::reset(double *vetor, int size){
	for(int i=0; i<size; i++){
		vetor[i] = 0;
	}
}

void variable::ls_nvalue(){
	
	find_nvalue();
	trail = allocate(nvalues);
	prob = allocate(nvalues);
}