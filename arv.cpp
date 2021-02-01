#include "arv.h"

arvore *criar_arvore(){
    return NULL;
}

arvore *criar_arvore(int valor, arvore *esq, arvore *dir){
    arvore *no = new arvore;
    no->valor = valor;
    no->esq = esq;
    no->dir = dir;
    return no;
}

void imprime_arvore(arvore *arv, int b){
    if(arv==NULL){
        for(int i=0; i<b; i++) cout << "\t";
        cout << "*\n";
        return;
    }
    imprime_arvore(arv->dir, b+2);
    imprimir(arv, b);
    imprime_arvore(arv->esq, b+2);
}

void imprimir(arvore *arv, int b){
    for(int i=0; i<b; i++) cout << "\t";
    cout << arv->valor << endl;
}

arvore *insert_partition(int inicio, int fim){
    int meio = (fim-inicio+1)/2+inicio;
    if(inicio==meio)
        return criar_arvore(inicio, NULL, NULL);
    arvore *arv = criar_arvore(meio, NULL, NULL);
    arv->esq = insert_partition(inicio, meio-1);
    if(meio < fim)
        arv->dir = insert_partition(meio+1, fim);
    return arv;
}

void liberar(arvore *arv){
    if(arv!=NULL){
        liberar(arv->esq);
        liberar(arv->dir);
        delete arv;
    }
}