// clases

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;


string string1 = "magnetizacion";
string string2 = "energia";
string string3 = "capacidad_calorifica";
string string4 = "suceptibilidad_magnética";


class variable_micro{
  
  
  private:
  
  vector<double> datos;
  
  public:
    

  variable_micro () {
  }
  
  void agregar_dato(double x){
    datos.push_back(x);
//     cout << datos[i] << end;
  }
  
  void reinicializar_vector(){
    datos.resize(0);
  }
  
  
 

  double promedio(){
    double suma=0;
//     cout << "numero de datos: " << datos.size() << endl;
    for (int i=0; i<(datos.size()); i++)  {
      suma+= datos[i];
    }
    return suma/(datos.size());
  }



  double promedio_cuadrado(){
    double suma=0;
    for (int i=0; i<(datos.size()); i++)  {
	suma+= datos[i]*datos[i];
      }
    return suma/(datos.size());
  }



  double varianza(){
    double sigma_cuadrada = (promedio_cuadrado()-(promedio())*(promedio()) ) ;
    return sigma_cuadrada;
  }
  
}; 

class variables_macro{
  private:
  
    vector<double> datos;
  
  public:

    variable_micro () {
    }
    
}; 
