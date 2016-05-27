#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <map>
#include "algoritmo.h"

using namespace std;


int main(int argc, char *argv[]){

    srand(time(0));
  
    double cambio_en_parametro;
    double campo_externo;
    double beta, T;
    char parametro;
    double f_min;
    int L, h, vecinos, frec_de_grafica, num_de_caminantes, repet_de_sim, tope_de_bloques;
    string red, modelo, algoritmo, longitud, caminantes, f_minimo, etiqueta_de_sim, etiqueta_de_bloques_max;
    map<string,int> valor_de_red;
    map<string,int> valor_de_modelo;
    map<string,int> valor_de_algoritmo;
    
    valor_de_red["cuadrada"]=1;
    valor_de_red["triangular"]=2;
    
    valor_de_modelo["Ising-ferro"]=1;
    valor_de_modelo["Ising-antiferro"]=2;
    
    valor_de_algoritmo["Metropolis"]=0;
    valor_de_algoritmo["WL1D"]=1;
    valor_de_algoritmo["WL2D"]=2;
    
    algoritmo = string(argv[3]);
    red = "_red-" + string(argv[1]);
    longitud = "_L-" + string(argv[4]);
    modelo = "_" +  string(argv[2]);
    caminantes = "_cam-" + string(argv[7]);
    f_minimo = "_f-" +  string(argv[6]);
    
    L = atoi(argv[4]);
    h = atoi(argv[5]);
    f_min = atof(argv[6]);
    num_de_caminantes = atoi(argv[7]);
    
    if (string(argv[8])=="no") {
        etiqueta_de_bloques_max = "";
    }
    else {
        etiqueta_de_bloques_max = "_max-bloques-" + string(argv[8]);
        tope_de_bloques = atoi(argv[8]);
    }
    
    if (string(argv[9])=="no") {
        frec_de_grafica = 0;
    }
    else {
        frec_de_grafica = atoi(argv[9]);
    }
    
    if (string(argv[10])=="no") {
        etiqueta_de_sim = "";
    }
    else {
        repet_de_sim = atoi(argv[10]);
        etiqueta_de_sim = "_sim-" + string(argv[10]);
    }
    
    string caracterizacion;
    caracterizacion = algoritmo + red + longitud + modelo + caminantes + f_minimo + etiqueta_de_bloques_max + etiqueta_de_sim;

    //
    cout << caracterizacion << endl;
    // Interpretando argumentos de la línea de comandos
    
    
    //Inicializando punteros a clases
    vector<Red*> r(num_de_caminantes);
    vector <Modelo*> m(num_de_caminantes);
    Algoritmo* a;
  
    // Definiendo punteros a red
    switch (valor_de_red[argv[1]]){
        case 1:
            
            for (int i=0; i<num_de_caminantes; i++) {
                r[i] = new Red_Cuadrada(L);
            }
            
            break;
        case 2:
            for (int i=0; i<num_de_caminantes; i++) {
                r[i] = new Red_Triangular(L);
            }
            break;
    }
    // extrayendo número de vecinos de cada espín
    vecinos = r[0]->numero_de_vecinos();
    // Definiendo punteros a modelo
    switch (valor_de_modelo[argv[2]]){
        case 1:
            for (int i=0; i<num_de_caminantes; i++) {
                m[i] = new Ising(r[i],h);
            }
            break;
        case 2:
            for (int i=0; i<num_de_caminantes; i++) {
                m[i] = new Antiferro(r[i],h,valor_de_red[argv[1]]);
            }
            break;
    }
    // Definiendo puntero a algoritmo
    switch (valor_de_algoritmo[argv[3]]){
        case 0:
            cout << "NADA" << endl;
            break;
        case 1:
            a = new Wang_Landau_1D(m, caracterizacion, num_de_caminantes, L, h, f_min, vecinos, frec_de_grafica);
            break;
        case 2:
            a = new Wang_Landau_2D(m, caracterizacion, num_de_caminantes, L, h, f_min, argv[8], vecinos,  frec_de_grafica);
            cout << "algoritmo creado\n";
            break;
    }
        //caracterizacion = algoritmo + red + "L_" + argv[4] + modelo + "_cam_" + argv[7] + "_f_" + argv[6];

}