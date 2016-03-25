#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "algoritmo.h"

using namespace std;


int main(int argc, char *argv[]){

    srand(time(0));
  
    
    
    double cambio_en_parametro;
    double HH;
    double beta, T;
    char parametro;
    int LL, vecinos;
    string modelo, red, algoritmo;
    
    
    // Interpretando argumentos de la línea de comandos
    int tipo_de_red;
    tipo_de_red = atoi(argv[1]);

    int tipo_de_modelo;
    tipo_de_modelo = atoi(argv[2]);
  
    int tipo_de_algoritmo;
    tipo_de_algoritmo = atoi(argv[3]);
    
  
    cout << "Dame el número de espines por lado: ";
    LL = atoi(argv[4]);
  
    cout << "Dame el valor del campo externo: ";
    HH = atoi(argv[5]);
    cout << HH << endl;
    
    double factor_de_mod_min;
    factor_de_mod_min = atof(argv[6]);
    cout << factor_de_mod_min << endl;
    
    int caminantes;
    caminantes = atoi(argv[7]);
    
    int animar;
    animar = atoi(argv[8]);
    
    int cantidad_de_sweeps;
    cantidad_de_sweeps = atoi(argv[9]);
    
    
    //Inicializando punteros a clases
    vector<Red*> r(caminantes);
    vector <Modelo*> m(caminantes);
    Algoritmo* a;
  
    // Definiendo punteros a red
    switch (tipo_de_red){
        case 1:
            red = "cuadrada_";
            for (int i=0; i<caminantes; i++) {
                r[i] = new Red_Cuadrada(LL);
            }
            
            break;
        case 2:
            red = "triangular_";
            for (int i=0; i<caminantes; i++) {
                r[i] = new Red_Triangular(LL);
            }
            break;
    }
  
    // extrayendo número de vecinos de cada espín
    vecinos = r[0]->numero_de_vecinos();
    
    
    // Definiendo punteros a modelo
    switch (tipo_de_modelo){
        case 1:
            modelo = "_Is_ferr";
            for (int i=0; i<caminantes; i++) {
                m[i] = new Ising(r[i],HH);
            }
            break;
        case 2:
            modelo = "_Is_antiferr";
            for (int i=0; i<caminantes; i++) {
                m[i] = new Antiferro(r[i],HH,tipo_de_red);
            }
            break;
    }
  
    string caracterizacion;
    // Definiendo puntero a algoritmo
    switch (tipo_de_algoritmo){
        case 1:
            cout << "NADA" << endl;
            break;
        case 2:
            algoritmo = "WL1D_";
            caracterizacion = algoritmo + red + "L_" + argv[4] + modelo + "_cam_" + argv[7] + "_f_" + argv[6];
            cout << caracterizacion << endl;
            a = new Wang_Landau_1D(m, caracterizacion, caminantes, LL, HH, factor_de_mod_min, vecinos, animar, cantidad_de_sweeps);
            break;
        case 3:
            algoritmo = "WL2D_";
            caracterizacion = algoritmo + red + "L_" + argv[4] + modelo + "_cam_" + argv[7] + "_f_" + argv[6];
            cout << caracterizacion << endl;
            a = new Wang_Landau_2D(m, caracterizacion, caminantes, LL, HH, factor_de_mod_min, vecinos, animar, cantidad_de_sweeps);
            break;
    }


}