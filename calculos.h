// obtención de datos

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "algoritmo.h"

using namespace std;


class Calcular{
  
  public:
    
    virtual void graficar()=0;
    //variables y funciones

    
};


class Resultados_W_L: public Calcular{
  
  private:

    ostringstream graf;
  
  public:
    Resultados_W_L(): {

      void graficar();

    }
    
    
    void graficar(){  
      graf.str("");
      graf << "set title 'Capacidad_calorifica'" << endl;
      graf << "plot \"Capacidad_calorifica_L_" << L << "_h_" << h << ".txt\" u 1:2 w lp" << endl;
//      graf << "set term wxt 3" << endl;
      graf << "set title 'Magnetizacion'" << endl;
      graf << "plot \"Magnetizacion_L_" << L << "_h_" << h << ".txt\" u 1:2 w lp" << endl;
//      graf << "set term wxt 4" << endl;
      graf << "set title 'Suceptibilidad_magnetica'" << endl;
      graf << "plot \"Suceptibilidad_magnetica_L_" << L << "_h_" << h << ".txt\" u 1:2 w lp" << endl;
//      graf << "set term wxt 5" << endl;
      graf << "set title 'Cumulante_de_binder'" << endl;
      graf << "plot \"Cumulante_de_binder_L_" << L << "_h_" << h << ".txt\" u 1:2 w lp" << endl;
//      graf << "set term wxt 6" << endl;
      graf << "set title 'Magnetizacion_staggered'" << endl;
      graf << "plot \"Magnetizacion_staggered_L_" << L << "_h_" << h << ".txt\" u 1:2 w lp" << endl;
//      graf << "set term wxt 7" << endl;
      graf << "set title 'Entropia_correcta'" << endl;
      graf << "plot \"Entropia_L_" << L << "_h_" << h << ".txt\" index 19 u 1:2 " << endl;
//      graf << "set term wxt 8" << endl;
      graf << "set title 'Entropia'" << endl;
      graf << "splot \"Entropia_WL1D_L_" << L << "_h_" << h << ".txt\" u 1:2:3 " << endl;
      FILE* gp;  // gnuplot es un puntero al objeto de tipo FILE
      gp = popen("gnuplot -persist", "w");
      fprintf(gp, "%s\n", graf.str().c_str());
      fflush(gp);
      cout << endl;
    }
    
};



