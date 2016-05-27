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

class Datos{

//   private:
    
    
  public:

    virtual void crear_archivos_de_datos()=0;

    
};


class Datos-Metropolis: public Datos{

//   private:

  public:

    virtual void crear_archivos_de_datos()=0;

};


class Datos_WL: public Datos{
  
  private:
    
    const int L;
    const int N;
    const int z; // número de vecinos
    const double h;
    
    
  public:  
    Red_Cuadrada(int LL, int H): L(LL), N(LL*LL), h(H){
        
        
    }
    
    void crear_archivos_de_datos(){
      nombre1.str("");
      nombre1 << "Entropia" << "_L_" << L << "_h_" << h << ".txt";
      archivo1.open(nombre1.str().c_str());
      nombre2.str("");
      nombre2 << "Histograma" << "_L_" << L << "_h_" << h << ".txt";
      archivo2.open(nombre2.str().c_str());
      nombre3.str("");
      nombre3 << "Capacidad_calorifica" << "_L_" << L << "_h_" << h << ".txt";
      archivo3.open(nombre3.str().c_str());
      nombre4.str("");
      nombre4 << "Magnetizacion" << "_L_" << L << "_h_" << h << ".txt";
      archivo4.open(nombre4.str().c_str());
      nombre5.str("");
      nombre5 << "Suceptibilidad_magnetica" << "_L_" << L << "_h_" << h << ".txt";
      archivo5.open(nombre5.str().c_str());
      nombre6.str("");
      nombre6 << "Cumulante_de_binder" << "_L_" << L << "_h_" << h << ".txt";
      archivo6.open(nombre6.str().c_str());
      nombre7.str("");
      nombre7 << "Magnetizacion_staggered" << "_L_" << L << "_h_" << h << ".txt";
      archivo7.open(nombre7.str().c_str());
      nombre8.str("");
      nombre8 << "Entropia_WL1D" << "_L_" << L << "_h_" << h << ".txt";
      archivo8.open(nombre8.str().c_str());
      nombre9.str("");
      nombre9 << "Entropia_1D" << "_L_" << L << "_h_" << h << ".txt";
      archivo9.open(nombre9.str().c_str());
    }
};
