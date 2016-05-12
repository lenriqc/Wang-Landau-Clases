#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "aleatorios.h"

using namespace std;

class Red{

//   private:
    
    
  public:

    virtual void geometria()=0;
    virtual void iniciar_red()=0;
    virtual void reordenar_red()=0;
    virtual void arreglo_de_subredes()=0;
    virtual int numero_de_espines_por_lado()=0;
    virtual int numero_de_espines()=0;

//    virtual void imprimir_arreglo()=0;
    virtual int numero_de_vecinos()=0;
    virtual int valor_del_espin(int a)=0;
    virtual int suma_de_vecinos(int a)=0;
    virtual void cambiar_valor_de_espin(int a)=0;
    virtual int valor_de_red_reordenada(int a)=0;
    
    
};


class Red_Cuadrada: public Red{
  
  private:
    vector<int> valores;
    vector<vector<int> > vecinos;
    vector<int> red_reordenada;
    vector<vector<int> > red_subred;
  
    
    const int L;
    const int N;
    const int z; // número de vecinos
    
    
  public:  
    Red_Cuadrada(int LL): L(LL), N(LL*LL), z(4){

      valores.resize(N);
      vecinos.resize(N);
      red_reordenada.resize(N);
      red_subred.resize(N/2);
      for (int i=0; i<N; i++){
	vecinos[i].resize(z);
      }
      for (int i=0; i<N/2; i++){
	red_subred[i].resize(2);
      }
      geometria();
      iniciar_red();
      
    }
    
    
      int numero_de_espines(){
	return N;
      }
      
      int numero_de_espines_por_lado(){
	return L;
      }
      
      int numero_de_vecinos(){
	return z;
      }
      
      
      void geometria(){
  
//	cout << "L = " << L << endl;
//	cout << "N = " << N << endl;
	int b;
	for(int i=0;i<N; i++){
	  b=(i+L)%L;
	  vecinos[i].resize(z);
	  
	  vecinos[i][0]=(b-1+L)%L+i-b;
	  vecinos[i][1]=(b+1+L)%L+i-b;
	  vecinos[i][2]=(i-L+2*N)%N;
	  vecinos[i][3]=(i+L)%N;
	}
      }
      
      
      void iniciar_red(){
	valores.resize(N);
	for (int i=0; i<N; i++){
	  double r=drand();
	  if (r<0.5)
	  {
	    valores[i]=1;
	  }
	  else
	  {
	    valores[i]=-1;
	  }
	}
      }

      
      void arreglo_de_subredes(){
	int k = -1;
	for (int i=0; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k; 
	  }
	}
	for (int i=0; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
      }
     

      void reordenar_red(){  
	int k = -1;
	for (int i=0; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k; 
	  }
	}
	for (int i=0; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
      }
      
      
      int valor_de_red_reordenada(int a){
	return red_reordenada[a];
      }

      
      void imprimir_arreglo(){
	for(int i=0; i<L; i++){
	  for(int j=0; j<L; j++){
	    if (valores[j+L*i]==1){
	      cout << "+ ";
	    }
	    else{
		cout << "- ";
	    }  
	  }
	  cout << endl;
	}
      }

  
      // regresa el valor del espin "a"
      int valor_del_espin(int a){
	return valores[a];
      }


      // suma los vecinos del espin "a"
      int suma_de_vecinos(int a){
	int suma=0;
	for (int j=0; j<(numero_de_vecinos()); j++){
	    suma += valores[vecinos[a][j]];
	}
	return suma;
      }
      
            
      // cambia valor del espin "a"
      void cambiar_valor_de_espin(int a){
	valores[a] = -valores[a];
      }

      
};




class Red_Triangular: public Red{
  
  private:
    
    vector<int> valores;
    vector<vector<int> > vecinos;  
    vector<int> red_reordenada;
    vector<vector<int> > red_subred;
    
    
    const int L;
    const int N;
    const int z; // número de vecinos
    
    
  public:  
    Red_Triangular(int LL): L(LL), N(LL*LL), z(6){

      valores.resize(N);
      vecinos.resize(N);
      red_reordenada.resize(N);
      red_subred.resize(N/3);
      for (int i=0; i<N; i++){
	vecinos[i].resize(z);
      }
      for (int i=0; i<N/3; i++){
	red_subred[i].resize(3);
      }
    
      geometria();
      iniciar_red();
      
    }
    
    
      int numero_de_espines(){
	return N;
      }
      
      int numero_de_espines_por_lado(){
	return L;
      }
      
      int numero_de_vecinos(){
	return z;
      }
      
      
      void geometria(){
  
//	cout << "L = " << L << endl;
//	cout << "N = " << N << endl;
	int b;
	int c;
	
	for(int i=0;i<N; i++){
	  b=(i+L)%L; // la coorenada x
	  c=(i)/L; // LA COODENADA Y
	  vecinos[i].resize(z);
	  
	  vecinos[i][0]=c*L + (b+1)%L;
	  vecinos[i][1]=c*L + (b-1+L)%L;
	  vecinos[i][2]=((c+1)%L)*L + b;
	  vecinos[i][3]=((c+1)%L)*L + (b+1)%L;
	  vecinos[i][4]=((c-1+L)%L)*L + b;
	  vecinos[i][5]=((c-1+L)%L)*L + (b-1+L)%L;
	}
      }
      
      
      void iniciar_red(){
	valores.resize(N);
	for (int i=0; i<N; i++){
	  double r=drand();
	  if (r<0.5)
	  {
	    valores[i]=-1;
	  }
	  else
	  {
	    valores[i]=1;
	  }
	}
      }
      
      void arreglo_de_subredes(){}
      
/*    void arreglo_de_subredes(){
	int k = -1;
	for (int i=0; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_subred[k]= L*i+j;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k; 
	  }
	}
	for (int i=0; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
      }*/
      
      void reordenar_red(){  //esto es la cuadrada reordenada, falta hacerla triangular
	int k = -1;
	for (int i=0; i<L;i+=2){
	  for (int j=0;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=2){
	  for (int j=2;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=2; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k; 
	  }
	}
	for (int i=0; i<L;i+=2){
	  for (int j=1;j<L;j+=2){
	    k++;
	    red_reordenada[L*i+j]= k; 
	  }
	}
	for (int i=1; i<L;i+=3){
	  for (int j=0;j<L;j+=3){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=2; i<L;i+=3){
	  for (int j=2;j<L;j+=3){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=0; i<L;i+=3){
	  for (int j=2;j<L;j+=3){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=1; i<L;i+=3){
	  for (int j=1;j<L;j+=3){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
	for (int i=2; i<L;i+=3){
	  for (int j=0;j<L;j+=3){
	    k++;
	    red_reordenada[L*i+j]= k;
	  }
	}
      }
      
      
      int valor_de_red_reordenada(int a){
	return red_reordenada[a];
      }
      
      
      
      void imprimir_arreglo(){
	for(int i=0; i<L; i++){
	  for(int j=0; j<L; j++){
	    if (valores[j+L*i]==1){
	      cout << "+ ";
	    }
	    else{
		cout << "- ";
	    }  
	  }
	  cout << endl;
	}
      }



      // regresa eñ valor del espin "a"
      int valor_del_espin(int a){
	return valores[a];
      }


      // suma los vecinos del espin "a"
      int suma_de_vecinos(int a){
	int suma=0;
	for (int j=0; j<(numero_de_vecinos()); j++){
	    suma += valores[vecinos[a][j]];
	}
	return suma;
      }
      
      
      
      // cambia valor del espin "a"
      void cambiar_valor_de_espin(int a){
	valores[a] = -valores[a];
      }



};
