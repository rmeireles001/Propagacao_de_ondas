#ifndef arv_h
#define arv_h

#include <iostream>

using namespace std;

struct arvore{
    int valor;
    arvore *esq, *dir;
};

arvore *criar_arvore();
arvore *criar_arvore(int valor, arvore *esq, arvore *dir);
void imprime_arvore(arvore *arv, int b);
void imprimir(arvore *arv, int b);
void imprimir(arvore *arv);
arvore *insert_partition(int inicio, int fim);
void liberar(arvore *arv);

#endif