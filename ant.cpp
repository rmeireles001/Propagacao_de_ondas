#include "ant.h"

ants::ants(){
	
}

void ants::init(int itno, double iofn){
	this->ofn = iofn;
	this->iter = itno;
	this->bestofn = 1e25;
}

void ants::copy_from(ants *from){
	ofn = from->ofn;
	bestofn = from->bestofn;
	iter = from->iter;
	variable = from->variable;
}