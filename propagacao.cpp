#include "propagacao.h"

propagacao::propagacao(){
	Gexp = new double[pontos];
	G = new double[pontos];
	A = new double[pontos];
	Z = new double[pontos];
	R = new double[pontos];
	F = new double[pontos];
	P = new double[pontos];
	config = new int[pontos];
	posicao = new double[pontos];
	A1 = alocar_mat(pontos);
	B1 = alocar_mat(pontos);
	Y = alocar_mat(pontos);
	for(int i=0; i<pontos; i++){
		Gexp[i] = 0.0;
		G[i] = 0.0;
		A[i] = 0.0;
		Z[i] = 0.0;
		R[i] = 0.0;
		F[i] = 0.0;
		P[i] = 0.0;
		posicao[i] = 0.0;
		config[i] = 0.0;
	}
	for(int i=0; i<pontos; i++){
		for(int j=0; j<pontos; j++){
			A1[i][j] = 0.0;
			B1[i][j] = 0.0;
			Y[i][j] = 0.0;
		}
	}
	E = 71*pow(10, 9);
	c = sqrt(E/ro);
	cout << "Alocados\n";
}

propagacao::~propagacao(){
	delete Gexp;
	delete G;
	delete A;
	delete Z;
	delete R;
	delete F;
	delete P;
	delete posicao;
	delete config;
	liberar_mat(A1, pontos);
	liberar_mat(B1, pontos);
	liberar_mat(Y, pontos);
	cout << "liberados\n";
}

double **propagacao::alocar_mat(int size){
	double **nv = new double*[size];
	for(int i=0; i<size; i++){
		nv[i] = new double[size];
	}
	return nv;
}

void propagacao::inserir(string s){
	int k;
	ifstream mf(s.c_str());
	for(int i=1; i<pontos; i++){
		mf >> posicao[i] >> A[i]; 
	}
	posicao[0] = 0;
	A[0] = A[1];
	/*for(int i=0, k=480; i<100; i++, k++){
		A[i] = A[k];
		posicao[i] = posicao[k];
	}*/
}

void propagacao::zerar(int n){
	for(int i=0; i<=n; i++){
		G[i] = 0.0;
		Z[i] = 0.0;
		R[i] = 0.0;
		F[i] = 0.0;
		P[i] = 0.0;
	}
	for(int i=0; i<=n; i++){
		for(int j=0; j<=n; j++){
			A1[i][j] = 0.0;
			B1[i][j] = 0.0;
			Y[i][j] = 0.0;
		}
	}
}

void propagacao::estimativa_inicial(double n){
	for(int i=0; i<=n; i++){
		A[i] = 300.0;
	}
}

void propagacao::config_area(double n){
	for(int i=0; i<=n; i++){
		config[i] = A[i]/A[0];
	}
}

void propagacao::attr_config(int k, double var){
	config[k] = var;
}

void propagacao::atualizar_area(int k){
	A[k] = config[k]*300;
}

void propagacao::dados_experimentais(string s, int n, int ini){
	int i, j;
	double *aux1, *aux2;
	double aux;
	ifstream mf(s.c_str());
	aux1 = new double[pontos];
	aux2 = new double[pontos];
	for(i=1; i<pontos; i++){
		mf >> aux1[i] >> aux >> aux2[i];
	}
	mf.close();
	for(i=ini, j=3; i<=(ini+n-3); i++, j++){
		Gexp[j] = aux2[i];
		posicao[j] = aux1[i];
	}
	Gexp[1] = aux2[ini-2];
	Gexp[2] = aux2[ini-1];
	posicao[1] = aux1[ini-2];
	posicao[2] = aux1[ini-1];
	delete aux1;
	delete aux2;
}

void propagacao::atribuirA(int n, double v){
	A[n] = v;
}

void propagacao::liberar_mat(double **mat, int size){
	for(int i=0; i<size; i++){
		delete mat[i];
	}
	delete mat;
}

void propagacao::prob_inverso(){
	int n1;
	F[1]=1.0;

	for(int i=0; i<=2; i++){
		Z[i]=300*A[i]*ro*c;
	}

	for(int i=1; i<=2; i++){
		if((Z[i]+Z[i-1])!=0.0)
			R[i] = (Z[i]-Z[i-1])/(Z[i]+Z[i-1]);
	}
	P[1] = 0.0;
	P[2] = 0.0;

	n1 = 1;
	G[1] = ((G[1] + (R[n1+1]+P[n1])*F[1-n1+1]));

	for(n1=1; n1<=2; n1++){
		G[2]=((G[2]+(R[n1+1]+P[n1])*F[2-n1+1]));
	}
}

void propagacao::prob_inverso(int n){
	double somatorio1;
	int n1, p1, k;

	Z[n]=300*A[n]*ro*c;

	if((Z[n]+Z[n-1])!=0.0)
		R[n] = (Z[n]-Z[n-1])/(Z[n]+Z[n-1]);
	for(p1=1; p1<=n-2; p1++){
		if(p1<=(n-3))
			A1[p1][n]=0;
		A1[n-2][n-1]=0;
		B1[p1][n-1] = 0;
		Y[p1][n] = 0;
		P[n] = 0;
	}

	for(p1=1; p1<=n-2; p1++){
		if(p1<=(n-3))
			A1[p1][n]=A1[p1][n-2]-B1[p1][n-2];
		A1[n-2][n-1]=0;
		somatorio1=0.0;
		for(k=1; k<=p1-1; k++)
			somatorio1=somatorio1+Y[k][n-1];
		B1[p1][n-1] = (R[n-p1-1]*(somatorio1+R[n-1]));
		Y[p1][n] = (R[n-p1]*(A1[p1][n-1] - B1[p1][n-1]));
		P[n] = (P[n]+Y[p1][n]);
	}

	G[n] = 0;
	for(n1=1; n1<=n; n1++){
		G[n]=((G[n]+(R[n1]+P[n1]))*F[n-n1+1]);
	}
}

void propagacao::prob_direto(int n, double *G){
	double somatorio1;
	int n1, p1, k;

	zerar(n);

	F[1]=1.0;

	for(int i=0; i<=n; i++){
		Z[i]=300*A[i]*ro*c;
	}

	for(int i=1; i<=n; i++){
		if((Z[i]+Z[i-1])!=0.0)
			R[i] = (Z[i]-Z[i-1])/(Z[i]+Z[i-1]);
	}

	P[1] = 0.0;
	P[2] = 0.0;

	n1 = 1;
	G[1] = ((G[1] + (R[n1+1]+P[n1])*F[1-n1+1]));

	for(n1=1; n1<=2; n1++){
		G[2]=((G[2]+(R[n1+1]+P[n1])*F[2-n1+1]));
	}

	for(n1=3; n1<=n; n1++){
		for(p1=1; p1<=n1-2; p1++){
			if(p1<=(n1-3))
				A1[p1][n1]=A1[p1][n1-2]-B1[p1][n1-2];
			A1[n1-2][n1-1]=0;
			somatorio1=0.0;
			for(k=1; k<=p1-1; k++)
				somatorio1=somatorio1+Y[k][n1-1];
			B1[p1][n1-1] = (R[n1-p1-1]*(somatorio1+R[n1-1]));
			Y[p1][n1] = (R[n1-p1]*(A1[p1][n1-1] - B1[p1][n1-1]));
			P[n1] = (P[n1]+Y[p1][n1]);
		}
	}
	for(int j=3; j<=n; j++){
		for(n1=1; n1<=j; n1++){
			G[j]=((G[j]+(R[n1]+P[n1]))*F[j-n1+1]);
		}
	}
}

double propagacao::erroG(int i){
	return sqrt(pow(G[i]-Gexp[i], 2));//pow(A[i]-0.5, 2)+2;//sqrt(pow(G[i]-Gexp[i], 2));
}

double prand(double low, double high)
{
	//srand(time(NULL));
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

double func(double x){
	return pow(x-2, 2)+2;
}

double propagacao::luus_jaakola(int pos){
	double mini, maxi, Rr, r, eps, aux1, aux2, aux1ant, qbest, t1, condicao, randomico, newconfig, oldconfig;
	int n_in, i, j, k, it, n_out, PAROU=0, aux, atrib, x, custo=0;
	Rr = 0.0;
	//G[] = 0.0;
	n_out = 1;
	n_in = 1;
	eps = 0.05;

	mini = 0.0;
	maxi = 1.0;

	r = maxi-mini;

	randomico = prand(0, 1);
	oldconfig = mini + r*randomico;

	aux1 = 1000;
	qbest = 1000;
	i = 1.0;
	condicao = pow(10, -10);
	while(qbest > condicao && PAROU<500){
		randomico = prand(0, 1);
		Rr = ((-0.5 + randomico)*(r));
		newconfig = oldconfig + Rr;
		if(newconfig < mini){
			newconfig = mini;
		}
		if(newconfig > maxi){
			newconfig = maxi;
		}
		atribuirA(pos, newconfig);
		prob_inverso(pos);
		aux1 = erroG(pos);
		custo++;
		atribuirA(pos, oldconfig);
		prob_inverso(pos);
		aux2 = erroG(pos);
		custo++;
		if(aux1<=aux2){
			qbest = aux1;
			oldconfig = newconfig;
		}

		i++;
		PAROU++;
		r = ((1 - eps)*r);
	}
	config[pos] = custo;
	return oldconfig;
}

void propagacao::run_lj(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	double img;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		img = luus_jaakola(i);
		atribuirA(i, img);
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "LUUS JAAKOLA", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
}

/*void propagacao::run_aco(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	aco a;
	a.get_data(0,1,0.05);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		a.run(i, this);
		atribuirA(i, a.get_var());
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "ANT COLONY OPTIMIZATION", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
}*/

//(double)(end-start)/(double)(CLOCKS_PER_SEC)
void propagacao::escrever_txt(string saida, string metodo, double tempo, int inicio, int fim){
	FILE *resultados = fopen(saida.c_str(), "w");
	fprintf(resultados, "Método usado: %s\n\nPOSIÇÃO\tÁREA(mm²)\tCUSTO\n\n", metodo.c_str());
	for(int i=inicio; i<=fim; i++){
		fprintf(resultados, "%d\t%.15e\t%d\n", (int) posicao[i], A[i], (int) config[i]);
	}
	fprintf(resultados, "\n\n");
	fprintf(resultados, "Tempo gasto: %lf\n", tempo);
	fclose(resultados);
}

double function(double x){
	return pow(x-1.5, 2)+3;
}

double propagacao::cgrasp(double llimit, double ulimit, double incr, int section){
	double fbest = 1e10;
	double xbest = 0;
	double f;
	//Busca linear / construção
	for(double x=llimit; x<=ulimit; x+=incr){
		atribuirA(section, x);
		prob_inverso(section);
		f = erroG(section);
		config[section]++;
		if(f<fbest){
			fbest = f;
			xbest = x;
		}
		cout << x << " " << ulimit << " " << ulimit+incr << endl;
	}
	//Busca local. Vizinhança anterior
	for(double x=xbest-incr; x>=llimit && x<xbest; x+=0.01){
		atribuirA(section, x);
		prob_inverso(section);
		f = erroG(section);
		config[section]++;
		if(f<fbest){
			fbest = f;
			xbest = x;
		}
	}
	//Busca local. Vizinhança posterior
	for(double x=xbest; x<=ulimit && x<=xbest+incr; x+=0.01){
		atribuirA(section, x);
		prob_inverso(section);
		f = erroG(section);
		config[section]++;
		if(f<fbest){
			fbest = f;
			xbest = x;
		}
	}
	return xbest;
}

void propagacao::run_cgrasp(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	double img;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		img = cgrasp(0, 1, 0.05, i);
		atribuirA(i, img);
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "C-GRASP", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
}
