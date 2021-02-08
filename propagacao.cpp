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
	condicao = pow(10, -5);
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
	double xi, xf;
	double f;
	double max=ulimit;
	int divisor = -2;
	//Busca linear / construção
	for(double x=llimit; x<=max && fbest>1e-5; x+=incr){
		atribuirA(section, x);
		prob_inverso(section);
		f = erroG(section);
		config[section]++;
		if(f<fbest){
			fbest = f;
			xbest = x;
		}
	}
	cout << endl;
	while(divisor>=-2){
		incr = pow(10, divisor);
		xi = xbest - incr*5;
		if(xi < llimit)
			xi = llimit;
		xf = xbest + incr*5;
		if(xf > ulimit)
			xf = ulimit;
		for(double x=xi; x<=xf && fbest>1e-5; x+=incr){
			atribuirA(section, x);
			prob_inverso(section);
			f = erroG(section);
			config[section]++;
			if(f<fbest){
				fbest = f;
				xbest = x;
			}
		}
		divisor--;
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

void propagacao::edit_areas(){
	for(int i=0; i<1000; i++){
		A[i] = 1.0;
		posicao[i] = (double)i+1;
	}
	/*double delta = 0.03;
	for(int i=150; i<165; i++){
		A[i] = A[i-1] - delta;
	}
	for(int i=165; i<175; i++){
		A[i] = A[i-1] + delta;
	}
	for(int i=175; i<180; i++){
		A[i] = A[i-1] + 0,05;
	}
	int j = 700;
	for(double k=-1; k<=-0.71; k+=0.01){
		j = ((int)(k*100))+800;
		A[j] = k*k;
		j++;
	}
	cout << j << endl;
	j = 729;
	for(double k=0; k<=0.2; k+=0.01){
		A[j] = 0.5+2.5*k;
		j++;
	}*/

	int j = 350;
	for(double k=0; k<=0.1; k+=0.01){
		A[j] = 1-5*k;
		j++;
	}
	for(double k=0; k<=0.1; k+=0.01){
		A[j] = 0.5+2.5*k;
		j++;
	}
	for(double k=0.06; k<=0.1; k+=0.01){
		A[j] = 0.5+5*k;
		j++;
	}
	cout << j << endl;
	j=650;
	for(double k=-0.28; k<=0.28; k+=0.01){
		A[j] = 0.5+6*k*k;
		j++;
	}
	cout << j << endl;
	FILE *f = fopen("areas2.txt", "w");
	for(int i=0; i<pontos-1; i++){
		fprintf(f, "%.15lf\t%.15e\n", posicao[i], A[i]);
	}
	fclose(f);
}

void propagacao::print_erro(int secao, double delta, string arq_ent, string arq_saida){
	FILE *saida = fopen(arq_saida.c_str(), "w");
	inserir(arq_ent);
	prob_direto(1000, Gexp);
	prob_direto(secao, G);
	fprintf(saida, "Posição: %d\n\nÁrea\tErro\n\n", secao);
	int max_value = 1/delta;
	/*for(double i=0; i<=1.001; i+=delta){
		atribuirA(secao, i);
		prob_inverso(secao);
		fprintf(saida, "%lf\t%lf\n", i, erroG(secao));
	}*/
	for(int sec=480; sec<=520; sec++){
		for(int i=0; i<=max_value; i++){
			atribuirA(sec, i*delta);
			prob_inverso(sec);
			fprintf(saida, "%lf\t", erroG(sec));
		}
		fprintf(saida, "\n");
	}
	fclose(saida);
}

double propagacao::binary_search(arvore *arv, double llimit, double ulimit, double delta, double secao){
	double g=0, gesq=10, gdir=10;
	double aresq=2, ardir=2, ar=2;
	int fesq=0, fdir=0;
	ar = arv->valor*delta+llimit;
	atribuirA(secao, ar);
	prob_inverso(secao);
	g = erroG(secao);
	if(arv->dir!=NULL){
		ardir = arv->dir->valor*delta+llimit;
		atribuirA(secao, ardir);
		prob_inverso(secao);
		gdir = erroG(secao);
		fdir = 1;
	}
	if(arv->esq!=NULL){
		aresq = arv->esq->valor*delta+llimit;
		atribuirA(secao, aresq);
		prob_inverso(secao);
		gesq = erroG(secao);
		fesq = 1;
	}
	/*cout << fesq << " " << (gesq<g) << " " << (aresq >= llimit) << endl;
	cout << fdir << " " << (gdir<g) << " " << (ardir<=ulimit) << endl;*/
	//cout << g << " " << gesq << " " << " " << gdir << endl;
	if(fesq && gesq<g && aresq>=llimit){
		if(fdir && gdir<gesq && ardir<=ulimit)
			return binary_search(arv->dir, llimit, ulimit, delta, secao);
		return binary_search(arv->esq, llimit, ulimit, delta, secao);
	}
	else if(fdir && gdir<g && ardir<=ulimit){
		return binary_search(arv->dir, llimit, ulimit, delta, secao);
	}
	return ar;
}

void propagacao::imprime_area(arvore *arv, int b, double delta, int secao, double llimit){
    if(arv==NULL){
        for(int i=0; i<b; i++) cout << "\t";
        cout << "*\n";
        return;
    }
    imprime_area(arv->dir, b+2, delta, secao, llimit);
    imprimir(arv->valor*delta+llimit, b, delta);
    imprime_area(arv->esq, b+2, delta, secao, llimit);
}

void propagacao::imprime_eco(arvore *arv, int b, double delta, int secao, double llimit){
    if(arv==NULL){
        for(int i=0; i<b; i++) cout << "\t";
        cout << "*\n";
        return;
    }
    imprime_eco(arv->dir, b+2, delta, secao, llimit);
	atribuirA(secao, arv->valor*delta+llimit);
	prob_inverso(secao);
    imprimir(erroG(secao), b, delta);
	imprimir(arv->valor*delta+llimit, b, delta);
    imprime_eco(arv->esq, b+2, delta, secao, llimit);
}

void propagacao::imprimir(double valor, int b, double delta){
    for(int i=0; i<b; i++) cout << "\t";
    cout << valor << endl;
}

double propagacao::busca_ordenada(ord *vet, int secao, double delta, double llimit, double ulimit){
	int max_value = (int)((ulimit-llimit)/delta);
	for(int i=0; i<=max_value; i++){
		vet[i].area = ((double) i)*delta+llimit;
		atribuirA(secao, vet[i].area);
		prob_inverso(secao);
		vet[i].ggexp = erroG(secao);
	}
	quickSort(vet, 0, max_value);
	return vet[0].area;
}

int propagacao::part(ord *vetor, int ini, int fim){
    int left, right;
	//double pivo;
	ord aux, pivo;
    left = ini;
    right = fim;
    pivo = vetor[ini];
    while(left<right){
        while(vetor[left].ggexp<=pivo.ggexp&&left<=fim){
            left++;
        }
        while(vetor[right].ggexp>pivo.ggexp&&right>=ini){
            right--;
        }
        if(left<right){
            aux = vetor[left];
            vetor[left] = vetor[right];
            vetor[right] = aux;
        }
    }
    vetor[ini] = vetor[right];
    vetor[right] = pivo;
    return right;
}

void propagacao::quickSort(ord *vetor, int ini, int fim){
    int pivo;
    if(fim>ini){
        pivo = part(vetor, ini, fim);
        quickSort(vetor, ini, pivo-1);
        quickSort(vetor, pivo+1, fim);
    }
}