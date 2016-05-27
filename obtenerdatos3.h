// obtención de datos

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "observables.h"
#include "algoritmo3.h"

using namespace std;


variable_micro magnetizacion;
variable_micro mstaggered;
variable_micro mtrianti;

variable_micro energia;

ostringstream nombre1, graf;

int N;
  double T, beta, h;  
  
void hacer_calculos(Algoritmo* a, Red* r, Modelo* m, char parametro, double dx){
  

  const int tiempo_total= 30000;
  const int numero_de_temperaturas= 40;
  const int tiempo_de_equilibracion = 10000;

  N = r->numero_de_espines();
  h = m->valor_de_h();
  nombre1.str("");
  nombre1 << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt";
  ofstream archivo1; // declarar "archivo" como un objeto de tipo ofstream // (output file)
  archivo1.open(nombre1.str().c_str());
  
//   if(parametro=='T'){
//       archivo1 << T << "\t\t";
//       archivo1 << "temp\t\t";
//     
//   }
//   else{
//     archivo1 << h << "\t\t";
//     archivo1  << "campo\t\t";
//   } 
//     
//   
//   archivo1 << "mag_stag\t\t" << "error de stag\t\t" << "mag alter\t\t" << "error de mag alter\t\t" << "magnetizacion\t\t" << "error de mag\t\t";
//   archivo1 << "energía\t\t" << "error en energía\t\t" << "capacidad calor\t\t" << "suceptibil magnet\t\t" << endl;
//     
  
  
  for (int j=0; j<numero_de_temperaturas; j++)  // este for me cambia la beta para crear tablas para distintas betas
  { 
    beta = a->valor_de_beta();
    T = 1/(beta);
    h = m->valor_de_h();
   
    if(parametro == 'T'){
      cout << T << endl;
    }
    else{
      cout << h << endl;
    }
  
    for (int i=0; i<tiempo_de_equilibracion; i++){
      a->sweep();
    }
    
    for (int i=tiempo_de_equilibracion; i<tiempo_total; i++){
      a->sweep();
      magnetizacion.agregar_dato(m->Magnetizacion());
      mstaggered.agregar_dato(m->MStaggered());
      mtrianti.agregar_dato(m->MTriAnti());
      
      energia.agregar_dato(m->Hamiltoniano());
    }
    
    
    
    archivo1.precision(5);
    
//     r->imprimir_arreglo();
    
    if(parametro=='T'){
      archivo1 << T << "\t\t";
//       archivo1 << "temp\t\t";
    
    }
    else{
      archivo1 << h << "\t\t";
//       archivo1  << "campo\t\t";
    } 
    
    archivo1 << magnetizacion.promedio() << "\t\t" << sqrt(magnetizacion.varianza()) << "\t\t";
    archivo1 << mstaggered.promedio() << "\t\t" << sqrt(mstaggered.varianza()) << "\t\t";
    archivo1 << mtrianti.promedio() << "\t\t" << sqrt(mtrianti.varianza()) << "\t\t";
    archivo1 << energia.promedio() << "\t\t" << sqrt(energia.varianza()) << "\t\t";
    archivo1 << (magnetizacion.varianza())*beta*N << "\t\t";
    archivo1 << (energia.varianza())*beta*beta*N << endl;
      
    magnetizacion.reinicializar_vector();
    energia.reinicializar_vector();
    mstaggered.reinicializar_vector();
    mtrianti.reinicializar_vector();
      
      
    if(parametro=='T'){
      a->aumentar_temperatura(dx);
    }
    else{
      m->aumentar_h(dx);
    }
  }
  
  archivo1.close();  
}



void grafica(){
  
    graf.str("");
    graf << "set title 'Magnetización tradicional'" << endl;
    graf << "plot \"" << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt" << "\" u 1:2\n";
    graf << "set term wxt 1\n";
    graf << "set title 'Magnetización Staggered'" << endl;
    graf << "plot \"" << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt" << "\" u 1:4 \n";
    graf << "set term wxt 2\n";
    graf << "set title 'Magnetización  propuesta'" << endl;
    graf << "plot \"" << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt" << "\" u 1:6\n";
    graf << "set term wxt 3\n";
    graf << "set title 'Energía'" << endl;
    graf << "plot \"" << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt" << "\" u 1:8\n";
    graf << "set term wxt 4\n";
    graf << "set title 'Suceptibilidad magnética'" << endl;
    graf << "plot \"" << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt" << "\" u 1:10\n";
    graf << "set term wxt 5\n";
    graf << "set title 'Capacidad calorífica'" << endl;
    graf << "plot \"" << "SIMULACION" << "_N_" << N << "_h_" << h << ".txt" << "\" u 1:11\n";
    FILE* gp;  // gnuplot es un puntero al objeto de tipo FILE
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "%s\n", graf.str().c_str());
    fflush(gp);
  }