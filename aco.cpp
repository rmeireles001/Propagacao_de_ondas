#include "aco.h"

void aco::get_data(double llimit, double ulimit, double increment){
	ant = new ants[ANTS];
	runbestant = new ants[runs];
	nevalbest = new long int[runs];
	runbestval = new values[runs];
	nevaluations=12;
	var.init(llimit, ulimit, increment);
	long int seed = time(NULL);
	#if localsearch
	ls = new lsaco;
	ls->init();
	#endif
}

void aco::ls_data(double llimit, double ulimit, double increment){
	ant = new ants[LSANTS];/*
	runbestant = new ants[runs];
	nevalbest = new long int[runs];
	runbestval = new values[runs];*/
	nevaluations=0;
	var.init(llimit, ulimit, increment);
	for(int antno=0; antno<LSANTS; antno++){
		ant[antno].init(0, 0.0);
	}
	bestant.init(0, 1e25);
	bestval.ofn=1e25;
}

void aco::end(){
	cout << "Destruindo geral" << endl;
	var.free();
	delete runbestant;
	delete nevalbest;
	delete runbestval;
	delete ant;
	#if localsearch
		ls->end();
		delete ls;
	#endif
}

void aco::initialize_ants_variables(int runno){
	nevaluations = 0;
	nevalbest[runno] = 0;
	for(int antno=0; antno<ANTS; antno++){
		ant[antno].init(0, 0.0);
	}
	bestant.init(0, 1e25);
	bestval.ofn=1e25;
}

void aco::iteration_init(int iterationno, int nants){
	for(int i=0; i<nants; i++){
		ant[i].iter = iterationno;
	}
	itbestant.init(iterationno, 1e25);
	itbestval.ofn = 1e25;
	itworstant.init(iterationno, 0.0);
}

//Trail
void aco::initialize_trail(){
	for(int j=0; j<var.nvalues; j++){
		var.trail[j] = tau0;
	}
}

void aco::evaporate_trail(){
	for(int j=0; j<var.nvalues; j++){
		var.trail[j] *= (1-rho);
	}
}

void aco::ras_update(double *help_b){
	int i, antno, target;
	double b;
	for(antno=0; antno<ANTS; antno++)
		help_b[antno] = ant[antno].ofn;
	for(i=0; i<ras_ranks-1; i++){
		b = help_b[0];
		target = 0;
		for(antno=0; antno<ANTS; antno++){
			if(help_b[antno]<b){
				b = help_b[antno];
				target = antno;
			}
		}
		help_b[target] = 1e25;
		update_trail_weighted(&ant[target], ras_ranks-1-i);
	}
	update_trail_weighted(&bestant, ras_ranks);
}

void aco::update_trail_weighted(ants *antptr, int weight){
	int vno, valueno;
	valueno=antptr->variable;
	if(antptr->ofn != 0.0){
		var.trail[valueno] += weight;
	}
	else{
		var.trail[valueno] += weight;
	}
}

void aco::trail(){
	double *help_b;
	evaporate_trail();
	t.reset(t.help_b, ANTS);
	ras_update(t.help_b);
}

void aco::find_values(int nants){
	for(int antno=0; antno<nants; antno++){
		ant[antno].variable = find_var_valueno();
	}
}

//find values for ants
int aco::find_var_valueno(){
	double q, sum, ranno;
	int valueno, nmax;
	find_prob();
	/*q = ran01(&seed);
	ranno = ran01(&seed);*/
	q = drand(0,1);
	ranno = drand(0,1);
	valueno = find_no_of_max(&nmax);

	if((q<q0) && (nmax==1)){
		#if(localupdate==1)
		var.trail[valueno] *= (1-gamma);
		#endif
		return valueno;
	}
	else{
		sum=0.0;
		for(valueno=0; valueno-1<var.nvalues && sum < ranno; valueno++){
			sum += var.prob[valueno];
		}
		#if(localupdate==1)
		if(valueno > 0 && valueno-1 < var.nvalues){
			var.trail[valueno-1] *= (1-rho);
		}
		#endif
		return (valueno-1);
	}
}

int aco::find_no_of_max(int *noofmax){
	int valueno, chooseno, wentinside=0, inloop=0;
	*noofmax=var.nvalues;
	double max=0;
	for(valueno=0; valueno<var.nvalues; valueno++){
		if(max < var.prob[valueno]){
			max = var.prob[valueno];
			chooseno = valueno;
		}
	}
	for(valueno=0; valueno<var.nvalues; valueno++){
		if(max > var.prob[valueno]){
			wentinside++;
			(*noofmax) = (*noofmax) - 1;
		}
	}
	return chooseno;
}

void aco::find_prob(){
	int valueno;
	double denominator=0.0;
	for(valueno=0; valueno<var.nvalues; valueno++){
		denominator += pow((double) var.trail[valueno], (double) alpha);
	}
	for(valueno=0; valueno<var.nvalues; valueno++){
		var.prob[valueno] = pow((double) var.trail[valueno], (double) alpha)/denominator;
	}
}

//Analysis
void aco::decode_varvalues(double *varvalues, int antno){
	int vno, valueno;
	if(ant[antno].variable == (var.nvalues)-1)
		*varvalues = var.ulimit;
	else
		*varvalues = (var.llimit+(ant[antno].variable*var.increment));
}

void aco::analysis(int iterationno, int nants, int section, propagacao *p){
	int antno, vno, valueno;
	t.tempval.iter = iterationno;
	for(antno=0; antno<nants; antno++){
		ant[antno].iter = iterationno;
		decode_varvalues(&t.tempval.var, antno);
		ant[antno].ofn = objective_function(&t.tempval.var, section, p);
		t.tempval.ofn = ant[antno].ofn;

		if(ant[antno].ofn < ant[antno].bestofn){
			ant[antno].bestofn = ant[antno].ofn;
			ant[antno].iter = iterationno;
		}

		if(ant[antno].ofn<itbestant.ofn){
			itbestval.copy_from(&t.tempval);
			itbestant.copy_from(&ant[antno]);
		}
		if(ant[antno].ofn>itworstant.ofn){
			itworstval.copy_from(&t.tempval);
			itworstant.copy_from(&ant[antno]);
		}
	}
}

double aco::objective_function(double *varvalue, int section, propagacao *p){
	p->atribuirA(section, *varvalue);
	p->prob_inverso(section);
	p->config[section]++;
	//p->prob_direto(section);
	//cout << "Varvalue: " << *varvalue << " Área: " << p->A[section] << " G: " << p->G[section] << " Gexp: " << p->Gexp[section] << " Erro: " << p->erroG(section) << endl;
	return p->erroG(section);
	//return pow(*varvalue-0.32, 2);
}

void aco::some_stats(int iterationno, int runno){
	if(itbestant.ofn<bestant.ofn){
		if(itbestval.ofn<bestval.ofn){
			bestant.copy_from(&itbestant);
			bestval.copy_from(&itbestval);
			nevalbest[runno] = nevaluations;
		}
	}
}

void aco::run(int section, propagacao *p){
	int runno, itno;
	double condicao = pow(10, -5);
	for(runno=0; runno<runs; runno++){
		initialize_ants_variables(runno);
		initialize_trail();
		for(itno=0; itno<ncmax && bestval.ofn > condicao; itno++){
			iteration_init(itno, ANTS);
			find_values(ANTS);
			analysis(itno, ANTS, section, p);
			some_stats(itno, runno);
			trail();
		}
		some_stats(itno, runno);
	}
	#if localsearch
	if(bestval.ofn > condicao){
		ls->run(&bestval, &var, section, p, &seed);
	}
	#endif
}

double aco::get_var(){
	return bestval.var;
}

void lsaco::init(){
	var.init(0, 0.1, 0.01);
	ant = (ants *) calloc(LSANTS, sizeof(ants));
	for(int antno=0; antno<LSANTS; antno++){
		ant[antno].init(0,0);
	}
}

void lsaco::get_data(values *gvalptr, variable *varv){
	double ll, ul, inc;
	ul = gvalptr->var + varv->increment;
	if(ul > varv->ulimit){
		ul = varv->ulimit;
	}
	ll = gvalptr->var - varv->increment;
	if(ll < varv->llimit){
		ll = varv->llimit;
	}
	var.increment = 0.01;
	var.ulimit = ul;
	var.llimit = ll;
	var.nvalues = 0;
	var.find_nvalue();
	for(int i=0; i<var.nvalues; i++){
		var.trail[i] = 0;
		var.prob[i] = 0;
	}
	for(int antno=0; antno<LSANTS; antno++){
		ant[antno].init(0,0);
	}
	bestant.init(0, 1e10);
}

void lsaco::run(values *gvalptr, variable *varv, int section, propagacao *p, long int *s){
	seed = *s;
	get_data(gvalptr, varv);
	initialize_trail();
	for(int itno=0; itno<lsncmax; itno++){
		iteration_init(itno, LSANTS);
		find_values(LSANTS);
		analysis(itno, LSANTS, section, p);
		some_stats(itno);
		trail();
	}
	if(bestant.ofn<gvalptr->ofn){
		gvalptr->copy_from(&bestval);
	}
}

void lsaco::some_stats(int iterationno){
	if(itbestant.ofn<bestant.ofn){
		bestant.copy_from(&itbestant);
		bestval.copy_from(&itbestval);
	}
}

void lsaco::evaporate_trail(){
	for(int j=0; j<var.nvalues; j++){
		var.trail[j] *= (1-lsrho);
	}
}

void lsaco::update_trail(ants *antptr){
	int vno, valueno;
	valueno=antptr->variable;
	if(antptr->ofn != 0.0){
		var.trail[valueno] += 1;
	}
	else{
		var.trail[valueno] += 1;
	}
}

void lsaco::trail(){
	evaporate_trail();
	for(int antno=0; antno<LSANTS; antno++){
		update_trail(&ant[antno]);
	//cout << "marcador" << endl;
	}
	update_trail(&itbestant);
}

void lsaco::end(){
	var.free();
	delete ant;
}

trash::trash(){
	help_b = new double[ANTS];
}

trash::~trash(){
	delete help_b;
}

void trash::reset(double *vec, int size){
	for(int i=0; i<size; i++){
		vec[i] = 0;
	}
}
double drand(double low, double high)
{
	//srand(time(NULL));
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

double ran01( long *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}