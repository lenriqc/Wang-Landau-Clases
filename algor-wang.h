#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "modelo.h"


using namespace std;

class Algoritmo{
  
  public:
    virtual void implementar_algoritmo()=0;
    virtual void sweep()=0;
    
    virtual double valor_de_beta()=0;
    virtual void aumentar_temperatura(double dx)=0;
    virtual void disminuir_temperatura(double dx)=0;
    
    virtual double calcular_magnetizacion()=0;
    virtual double calcular_energia()=0;


};

class Metropolis-Wang: public Algoritmo{
  
  private:
    Modelo* M;
    int N;
    vector <double> conjunto_de_energias;
    vector <double> densidad_de_estados;
    vector <int> histograma;
    
  public:
    
    Metropolis-Wang (Modelo* mod){
      M = mod;
      N = M->numero_de_espines();
      conjunto_de_energias.resize(N);
      densidad_de_estados.resize(N);
      histograma.resize(N);
    }
    
  void densidad_de_estados_inicial(){
  
    for (int i=0; i<N; i++){
      densidad_de_estados[i]=1;
    }
  }
  
  double densidad_de_estados(double E){
    
  
  }
  
  void obtener_energias_posibles(){
    conjunto_de_energias[0]= M->Hamiltoniano;
//     for (int i=1; i<N; i++){
  }
  
  void implementar_algoritmo(){
    int a = M->elegir_espin();
    double dH = M->cambio_en_energia(a);
    double Energia_inicial, Energia_final;
    Energia_inicial = M->Hamiltoniano;
    Energia_final = Energia_inicial + dH;
    double p;
    if (dH <= 0) {
      p=1;
      M->aceptar_cambio(a);
    }
    else {
      p = exp(-beta*dH);
      double r = drand();
      if (r<=p){
        M->aceptar_cambio(a);
      }
    }
  }
      
      
      
      
      
      
      
      