#include "values.h"

values::values(){
	var = 0;
}

void values::copy_from(values *from){
	ofn = from->ofn;
	iter = from->iter;
	var = from->var;
}