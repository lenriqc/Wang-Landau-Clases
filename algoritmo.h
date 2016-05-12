#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "modelo.h"


//Revisar
// 1 en estimar_entropia_2D la escritura en los archivos va después de renormalizar OJO!
// En la función que escribe a archivo final, se tiene que usar entropia para verificar la escritura, ya que si se usa el histograma nunca escribirá nada ya que el histograma se resetea.
// Tal vez no es necesario definir contadores_1D, parametro_1D, etc en el algor WL-1D, hay que preguntar. Lo mismo pero para 2D pasaría en algor WL-2D
//

using namespace std;

class Algoritmo{
    
  
    public:
    

//    virtual void implementar_algoritmo()=0;
//    virtual void sweep()=0;
    
/*    virtual int calcular_energia_configuracional()=0;
    virtual int calcular_magnetizacion()=0;
    virtual int calcular_parametro_de_orden()=0;*/

    virtual void crear_archivos_de_datos()=0;
//    virtual void datos_de_histograma()=0;
    virtual void obtener_variables_termodinamicas(int a)=0;
//    virtual void imprimir_resultados()=0;
    virtual void graficar()=0;

    };



class Metropolis: public Algoritmo{
  
  private:
    Modelo* M;
    double beta, T;
    int N;
    
  public:

    Metropolis(Modelo* mod, double Beta): beta(Beta) {
      M = mod;
      N = M->numero_de_espines();
      beta=1;
      T = 1/beta;
    }



    void implementar_algoritmo(){
      int a = M->elegir_espin();
      double dE = M->cambio_en_energia_configuracional(a);
      double p;
      if (dE <= 0) {
	p=1;
	M->aceptar_cambio(a);
      }
      else {
	p = exp(-beta*dE);
	double r = drand();
	if (r<=p){
	  M->aceptar_cambio(a);
	}
      }
    }  
        
    void sweep()
    {
      for (int i=0; i<N; i++)
      {
	implementar_algoritmo();
      }
    }
    
    
    void aumentar_temperatura(double delta_parametro){
      T += delta_parametro;
      beta = 1/T;
    }
    
    void disminuir_temperatura(double delta_parametro){
      T -= delta_parametro;
      beta = 1/T;
    }
    
    double valor_de_beta(){
      return(beta);
    }
    
    
    double calcular_magnetizacion(){
      return (M->magnetizacion());
    }
    
    
    double calcular_energia_configuracional(){ //checar cambio de nombre en función
      return (M->energia_configuracional()); //ya no hay hamiltoniano
    }
    
    
    double calcular_parametro_de_orden(){
      return (M->parametro_de_orden());
    }
    
    
    void crear_archivos_de_datos(){
    }
    
    void obtener_variables_termodinamicas(double beta){
    }
    
    void imprimir_resultados(){
    }

  
};


class Wang_Landau_1D: public Algoritmo{
    private:
    int caminantes;
    vector<Modelo*> M;
    string especificacion;
    int N;
    int L;
    int h;
    int t;
    int v;
    int E_min;
    int s;
    int paso;
    double f, f_minima;
    int a;
    int J1, J2;
    double numero_de_refinamientos;
    double logaritmo_de_factor_de_modificacion;
    double beta;
    double temperatura;
    double constante_de_renormalizacion;
    int cantidad_de_datos;
    int anima;
    int numero_de_sweeps;
    int contador_de_sweeps;
    int frecuencia_de_animacion;
    
    int estatus;
    string nombre[7];
    ofstream archivo[7];
    ostringstream flujo[7];
 
    FILE* gp[2];  // gnuplot es un puntero al objeto de tipo FILE
    vector <vector<double> > entropia_2D;
    vector <double> entropia_1D;
    
    vector <int> histograma_1D;
    vector <vector<int> > histograma_2D;
    vector <vector<double> > probabilidad_condicionada;
    
    vector <vector<int> > contador_de_parametro_de_orden;
    vector <vector<double> > parametro_de_orden_microcanonico;
    vector <vector<double> > parametro_de_orden_microcanonico_promedio;
    vector <int> contador_de_parametro_de_orden_1D;
    vector <double> parametro_de_orden_microcanonico_1D;
    vector <double> parametro_de_orden_microcanonico_promedio_1D;

    
    vector <vector<double> > resultados;
    
    ostringstream graf, anim;
    
    
    public:
    Wang_Landau_1D(vector<Modelo*> mod, string descripcion, int Caminantes, int l_de_red, double h_exter,  double factor_min, int vecinos, int frec_de_anim): h(h_exter ) {
        especificacion = descripcion;
        f_minima = factor_min;
        frecuencia_de_animacion = frec_de_anim;
        numero_de_sweeps = 1000;
        contador_de_sweeps = 0;
        caminantes = Caminantes;
        M.resize(caminantes);
        for (int i=0; i<caminantes; i++) {
            M[i] = mod[i];
        }
        N = M[0]->numero_de_espines();
        logaritmo_de_factor_de_modificacion = 1;
        L = l_de_red;
        N = L*L;
        f = 1;
        numero_de_refinamientos = 0;
        temperatura = 0.5;
        cantidad_de_datos = 50;
        J1 = M[0]->valor_de_J1();
        J2 = M[0]->valor_de_J2();
        paso=4;
        v = vecinos/2; //no está en uso
        E_min = (vecinos/2)*N;
        s = E_min/4;
        t = vecinos*N/paso + 1;
        if (frecuencia_de_animacion != 0) {
            anima = 1;
        }
        else {
            anima = 0;
        }
        //VECTORES
        entropia_1D.resize(t);
        histograma_1D.resize(t); //no aparece en wl2d
        entropia_2D.resize(t);
        histograma_2D.resize(t);
        probabilidad_condicionada.resize(t);// no aparece en wl2d
        contador_de_parametro_de_orden.resize(t);
        parametro_de_orden_microcanonico.resize(t);
        parametro_de_orden_microcanonico_promedio.resize(t);
        contador_de_parametro_de_orden_1D.resize(t);
        parametro_de_orden_microcanonico_1D.resize(t);
        parametro_de_orden_microcanonico_promedio_1D.resize(t);
        //NOMBRES
        nombre[0]="Histograma_animado";
        nombre[1]="Histograma_1D";
        nombre[2]="Histograma_2D";
        nombre[3]="Entropia_1D";
        nombre[4]="Entropia_2D";
        nombre[5]="Entropia_1D_final";
        nombre[6]="Entropia_2D_final";
        //CONSTRUCCIÓN DE VECTORES 2D
        for (int i=0;i<t;i++){
            histograma_2D[i].resize(N+1);
            entropia_2D[i].resize(N+1);
            probabilidad_condicionada[i].resize(N+1); // no está en wl2d
            contador_de_parametro_de_orden[i].resize(N+1);
            parametro_de_orden_microcanonico[i].resize(N+1);
            parametro_de_orden_microcanonico_promedio[i].resize(N+1);
        }
        //
        resultados.resize(5);
        for (int i=0;i<5;i++){
            resultados[i].resize(cantidad_de_datos);
        }
        //INICIALIZANDO VECTORES
        reiniciar_entropia_1D();
        reiniciar_histograma_1D(); //no aparece en WL2D
        reiniciar_entropia_2D();
        reiniciar_histograma_2D();
        reiniciar_parametro_de_orden_microcanonico();
        reiniciar_parametro_de_orden_microcanonico_promedio();
        reiniciar_contador_de_parametro_de_orden();
        reiniciar_parametro_de_orden_microcanonico_1D();
        reiniciar_parametro_de_orden_microcanonico_promedio_1D();
        reiniciar_contador_de_parametro_de_orden_1D();
        reiniciar_probabilidad_condicionada(); // no aparece en WL2D
        //ABRIR PIPE A GNUPLOT
        iniciar_conexion_a_gnuplot();
        //CREAR Y ABRIR ARCHIVOS PARA GUARDAR RESULTADOS
        crear_archivos_de_datos();
        //DATOS INICIALES DE CAMINANTES
        for (int i=0; i<caminantes; i++) {
            M[i]->Energia = M[i]->energia_configuracional();
            M[i]->Magnetizacion = M[i]->magnetizacion();
            M[i]->Parametro_de_orden = M[i]->parametro_de_orden();
        }
        for (int i=0; i<caminantes; i++) {
            cout << "E" << i << ": " <<M[i]->Energia << "\t" << "M" << i << ": " << M[i]->Magnetizacion << "\t" << "PO" << i << ": " << M[i]->Parametro_de_orden << endl;
        }
        //
        estimar_entropia_1D();//NO ESTA EN WL2D Y FALTA LA PROB COND
        //
        calcular_parametro_de_orden_microcanonico_promedio_1D();
        //
        calcular_parametro_de_orden_microcanonico_promedio();
        //
        calcular_probabilidad_condicionada();
        //
        estimar_entropia_2D();
        //
        entropia_1D_a_archivo();
        //
        entropia_2D_a_archivo();
        //
        //imprimir_resultados();
        //
        graficar();
        
        remove(flujo[1].str().c_str());
        remove(flujo[2].str().c_str());
        remove(flujo[3].str().c_str());
        remove(flujo[4].str().c_str());
    }
    
    
    void crear_archivos_de_datos(){
        for (int i=0; i<7; i++) {
            flujo[i].str("");
            flujo[i] << nombre[i] << "/" << nombre[i] << "_" << especificacion << ".txt";
            estatus=mkdir(nombre[i].c_str(),S_IRWXU);
            if (!estatus) {
                cout << "Carpeta " << nombre[i] << " creada." << endl;
            }
            else{
                cout << "Carpeta " << nombre[i] << " ya existe." << endl;
            }
            archivo[i].open(flujo[i].str().c_str());
        }
    }
    
    
    void iniciar_conexion_a_gnuplot(){
        anim.str("");
//        anim << "set term wxt 0" << endl;
        anim << "set xrange [" << -E_min << ":" << E_min << "]" << endl;
        anim << "set title 'Histograma animado'" << endl;
        gp[0] = popen("gnuplot -persist", "w");
        //
        graf.str("");
//        graf << "set term wxt 1 persist" << endl;
        gp[1] = popen("gnuplot -persist", "w");
    }
    
    
    void grafica_evolucion_de_histograma(){
        anim << "plot \"" << flujo[0].str().c_str() << "\"" << endl;
        fprintf(gp[0], "%s\n", anim.str().c_str());
        anim.str("");
        fflush(gp[0]);
    }

    
    void datos_de_histograma(){
        archivo[0].open(flujo[0].str().c_str());
        archivo[0] << "#" <<endl;
        for (int i=0;i<t;i++){
            if (histograma_1D[i]!=-1) {
                archivo[0] << 4*i-E_min << "\t" << histograma_1D[i] << endl;
            }
        }
        archivo[0].close();
    }

    
    void reiniciar_entropia_1D(){
        for (int i=0;i<t;i++){
            entropia_1D[i] = 0;
        }
    }  


    void reiniciar_entropia_2D(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                entropia_2D[i][j] = 0;
            }
        }
    }


    void reiniciar_histograma_1D(){
        for (int i=0;i<t;i++){
            histograma_1D[i] = -1;
        }
    }


    void reiniciar_histograma_2D(){
      for (int i=0;i<t;i++){
          for (int j=0;j<N+1;j++){
              histograma_2D[i][j] = -1;
          }
      }
    }


    void reiniciar_probabilidad_condicionada(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                probabilidad_condicionada[i][j] = -1;
            }
        }
    }


    void reiniciar_contador_de_parametro_de_orden(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                contador_de_parametro_de_orden[i][j] = -1;
            }
        }
    }


    void reiniciar_parametro_de_orden_microcanonico(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                parametro_de_orden_microcanonico[i][j] = -1;
            }
        }
    }


    void reiniciar_parametro_de_orden_microcanonico_promedio(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                parametro_de_orden_microcanonico_promedio[i][j] = -1;
            }
        }
    }
    
    
    void reiniciar_contador_de_parametro_de_orden_1D(){
        for (int i=0;i<t;i++){
            contador_de_parametro_de_orden_1D[i]= -1;
        }
    }
    
    
    void reiniciar_parametro_de_orden_microcanonico_1D(){
        for (int i=0;i<t;i++){
            parametro_de_orden_microcanonico_1D[i]= -1;
        }
    }
    
    
    void reiniciar_parametro_de_orden_microcanonico_promedio_1D(){
        for (int i=0;i<t;i++){
            parametro_de_orden_microcanonico_promedio_1D[i] = -1;
        }
    }


    void alterar_logaritmo_de_factor_de_modificacion(){
        f = f/2;
    }
    
    
    void proponer_cambio_configuracional(Modelo *MOD){
          a = MOD->elegir_espin();
          MOD->Energia_nueva = MOD->Energia + MOD->cambio_en_energia_configuracional(a);
          MOD->Magnetizacion_nueva = MOD->Magnetizacion + MOD->cambio_en_magnetizacion(a);
          MOD->Parametro_de_orden_nuevo = MOD->Parametro_de_orden + MOD->cambio_en_parametro_de_orden(a);
    }
    
    
    void implementar_algoritmo(Modelo *MOD){
        double dS = entropia_1D[(MOD->Energia_nueva)/4 + s] - entropia_1D[(MOD->Energia)/4 + s];
        double p;
        if (dS <= 0) {
            p=1;
            MOD->aceptar_cambio(a);
            MOD->actualizar_variables();
        }
        else {
            p = exp(-dS);
            double r = drand();
            if (r<=p){
                MOD->aceptar_cambio(a);
                MOD->actualizar_variables();
            }
            else {
            }
        }
        int k = (MOD->Energia)/4 + s;
        int l = (MOD->Magnetizacion)/2 + N/2;
        entropia_1D[k] += f;
        histograma_1D[k]++;
        contador_de_parametro_de_orden_1D[k]++;
        parametro_de_orden_microcanonico_1D[k] += abs(MOD->Parametro_de_orden);
        contador_de_parametro_de_orden[k][l]++;
        parametro_de_orden_microcanonico[k][l] += abs(MOD->Parametro_de_orden);
    }
  

    void sweep(){
        for (int i=0;i<N/caminantes;i++){
            for (int j=0; j<caminantes; j++) {
                proponer_cambio_configuracional(M[j]);
                implementar_algoritmo(M[j]);
            }
        }
    }

    
    void varios_sweeps(){
        for (int i=0; i<numero_de_sweeps; i++) {
            sweep();
            contador_de_sweeps++;
            if (contador_de_sweeps<frecuencia_de_animacion) {
            }
            else {
                contador_de_sweeps = 1;
                if (anima==1) {
                    datos_de_histograma();
                    grafica_evolucion_de_histograma();
                }
            }
        }
    }

// El último i=t siempre será cero y no se renormaliza bien =(
    void renormalizar_entropia_1D(){
        double min = 0;
        int k=0;
        while (entropia_1D[k] == 0) {
            k++;
        }
        min = entropia_1D[k];
        for (int i=k;i<t;i++){
            if (entropia_1D[i] != 0) {
                if (entropia_1D[i] <= min) {
                    min = entropia_1D[i];
                }
            }
        }
        for (int i=0;i<t;i++){
            if (entropia_1D[i] != 0) {
                entropia_1D[i] -= min-1;//para que el último i no se anule se pone el -1
            }
            
        }
    }


    void estimar_entropia_1D(){
        while (f>f_minima){
            varios_sweeps();
            cout << "f: " << f << endl;
            int k=0;
            while (k<t){
                if ((histograma_1D[k] >= 10) || (histograma_1D[k] == -1)) {
                    k++;
                }
                else {
                    varios_sweeps();
                }
            }
            renormalizar_entropia_1D();
            archivo[3] << endl << endl;
            archivo[1] << endl << endl;
            for (int i=0; i<t; i++){
                if (entropia_1D[i] !=0) {
                    archivo[3] << 4*i-E_min << "\t" << entropia_1D[i] << endl;
                    archivo[1] << 4*i-E_min << "\t" << histograma_1D[i] << endl;
                }
            }
            reiniciar_histograma_1D();
            alterar_logaritmo_de_factor_de_modificacion();
            numero_de_refinamientos++;
        }
    }
    
    
//REVISAR SI SE TIENE QUE USAR MAX O NORMA
    void calcular_probabilidad_condicionada(){
        double norma=0;
        for (int i=0;i<t;i++){
//            cout << "Para E = " << 4*i-E_min << endl;
            for (int j=0;j<N+1;j++){
                if (contador_de_parametro_de_orden[i][j] != -1){
                    norma += contador_de_parametro_de_orden[i][j]+1;//aquí había un +1 que al parecer no va, creo que sí va
                }
            }
            for (int j=0;j<N+1;j++){
                if (contador_de_parametro_de_orden[i][j] != -1){
                    probabilidad_condicionada[i][j] = contador_de_parametro_de_orden[i][j] / norma;
                }
            }
            norma = 0;
        }
    }
    

    void estimar_entropia_2D(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if (probabilidad_condicionada[i][j] != -1){
                    entropia_2D[i][j] = entropia_1D[i] + log(probabilidad_condicionada[i][j]);
                    archivo[4] << 4*i-E_min << "\t" << 2*j-N << "\t" << entropia_2D[i][j] << endl;
                }
            }
        }
    }
    
    
    void calcular_parametro_de_orden_microcanonico_promedio_1D(){
        for (int i=0;i<t;i++){
            if ( contador_de_parametro_de_orden_1D[i] != -1){
                parametro_de_orden_microcanonico_promedio_1D[i] = ((parametro_de_orden_microcanonico_1D[i] + 1.)/(contador_de_parametro_de_orden_1D[i] + 1.));
                cout << "PO("<< 4*i-E_min << "): " << parametro_de_orden_microcanonico_promedio_1D[i] << endl;
            }
        }
    }
    
    
    void calcular_parametro_de_orden_microcanonico_promedio(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if ( contador_de_parametro_de_orden[i][j] != -1){
                    parametro_de_orden_microcanonico_promedio[i][j] = ((parametro_de_orden_microcanonico[i][j] + 1.)/(contador_de_parametro_de_orden[i][j] + 1.));
                }
            }
        }
    }
    
    
    void entropia_1D_a_archivo(){
        archivo[5] << "Energía \t Magnetización \t Entropía \t P.O." << endl;
        for (int i=0;i<t;i++){
            if (entropia_1D[i]!=0) {
                archivo[5] << 4*i-E_min << "\t" << entropia_1D[i] << "\t" <<parametro_de_orden_microcanonico_promedio_1D[i] << endl;
            }
        }
        archivo[5].close();
    }
    
    
    void entropia_2D_a_archivo(){
        archivo[6] << "Energía \t Magnetización \t Entropía \t P.O." << endl;
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if ( entropia_2D[i][j] != 0){
                    archivo[6] << 4*i-E_min << "\t" << 2*j-N << "\t" << entropia_2D[i][j] << "\t" << parametro_de_orden_microcanonico_promedio[i][j] << endl;
                }
            }
        }
        archivo[6].close();
    }

    
    void encontrar_constante_de_renormalizacion(){
        constante_de_renormalizacion = 0;
        int E=0;
        int M=0;
        for (int i=0;i<t;i++){
            E = 4*i-E_min;
            for (int j=0;j<N+1;j++){
                if ((entropia_2D[i][j] - beta*(J1*E))>=constante_de_renormalizacion){
                    M = 2*j-N;
                    constante_de_renormalizacion = entropia_2D[i][j] - beta*(J1*E-J2*h*M);
                }
            }
        }
    }
    
    
    void obtener_variables_termodinamicas(int k){
        double funcion_de_particion = 0;
        double funcion_de_particion_2D = 0;
        double energia_promedio = 0;
        double energia_promedio_cuadrada = 0;
        double magnetizacion_promedio = 0;
        double magnetizacion_cuadrada_promedio = 0;
        double magnetizacion_cuarta_promedio = 0;
        double parametro_de_orden_promedio = 0;
        double exponente=0;
        int E=0;
        int M=0;
        for (int i=0;i<t;i++){
            if(entropia_1D[i] != 0){
                E = 4*i-E_min;
                for (int j=0;j<N+1;j++){
                    if(probabilidad_condicionada[i][j] != -1){
                        M = 2*j-N;
                        exponente = entropia_2D[i][j] - (beta)*(E-J2*h*M) - constante_de_renormalizacion;
                        funcion_de_particion_2D += exp(exponente);
                        energia_promedio += (abs(E))*exp(exponente);
                        energia_promedio_cuadrada += (abs(E*E))*exp(exponente);
                        magnetizacion_promedio += (abs(M))*exp(exponente);
                        magnetizacion_cuadrada_promedio += (abs(M*M))*exp(exponente);
                        magnetizacion_cuarta_promedio += pow(M,4.)*exp(exponente);
                        parametro_de_orden_promedio += (parametro_de_orden_microcanonico_promedio[i][j])*exp(exponente);
                    }
                }
            }
        }
        resultados[0][k] = (energia_promedio_cuadrada/funcion_de_particion_2D - (energia_promedio*energia_promedio)/(funcion_de_particion_2D*funcion_de_particion_2D))*(1./N) ;
        resultados[1][k] = magnetizacion_promedio/(N*funcion_de_particion_2D);
        resultados[2][k] = (magnetizacion_cuadrada_promedio/funcion_de_particion_2D - (magnetizacion_promedio*magnetizacion_promedio)/(funcion_de_particion_2D*funcion_de_particion_2D))/N;
        resultados[3][k] = 1 - ((magnetizacion_cuarta_promedio*funcion_de_particion_2D)/(3*magnetizacion_cuadrada_promedio*magnetizacion_cuadrada_promedio));
        resultados[4][k] = parametro_de_orden_promedio/(N*funcion_de_particion_2D);
        
    }
    
    
/*    void imprimir_resultados(){
        arch_cap_cal << endl << endl;
        arch_magnet << endl << endl;
        arch_sucep_magnet << endl << endl;
        arch_cbinder << endl << endl;
        arch_po << endl << endl;
        for (int i=0;i<cantidad_de_datos;i++){
            beta = 1/temperatura;
            //encontrar_constante_de_renormalizacion();
            obtener_variables_termodinamicas(i);
            arch_cap_cal << temperatura << "\t" << resultados[0][i];
            arch_cap_cal << endl;
            arch_magnet << temperatura << "\t" << resultados[1][i];
            arch_magnet << endl;
            arch_sucep_magnet << temperatura << "\t" << resultados[2][i];
            arch_sucep_magnet << endl;
            arch_cbinder << temperatura << "\t" << resultados[3][i];
            arch_cbinder << endl;
            arch_po << temperatura << "\t" << resultados[4][i];
            arch_po << endl;
            temperatura += 0.1;
        }
    //  cout << "número de refinamientos: " << numero_de_refinamientos -1 << endl;
    }*/
    
    void graficar(){
        graf << "encabezado = system('head -1 " << flujo[6].str().c_str() << " ')" << endl;
        graf << "set xlabel word(encabezado,1)" << endl;
        graf << "set ylabel word(encabezado,2)" << endl;
        graf << "set zlabel word(encabezado,3)" << endl;
        graf << "splot \"" << flujo[6].str().c_str() << "\" using 1:2:3" << endl;
        fprintf(gp[1], "%s\n", graf.str().c_str());
        fflush(gp[1]);
/*        cout << endl;
        int espera;
        cin >> espera;*/
    }

};





class Wang_Landau_2D: public Algoritmo{
    private:
    int caminantes;
    vector<Modelo*> M;
    string especificacion;
    int N;
    int L;
    int h;
    int paso; //paso mínimo entre dos energías
    int t; //tamaño del vector
    int s;
    int v;
    double f, f_minima;
    int a;
    int E_min;
    int J1, J2;
    double numero_de_refinamientos, logaritmo_de_factor_de_modificacion;
    double beta, temperatura;
    double constante_de_renormalizacion;
    int cantidad_de_datos;
    int tiempo_de_simulacion, contador_de_sweeps, contador_de_bloques, num_max_de_bloques, contador_de_animacion;
    string cantidad_max_de_bloques;
    int anima;
    int frecuencia_de_animacion;
    int condicion_de_histograma;
    bool proposicion_de_bloques_max;
    FILE* gp[2];
    int estatus;
    string nombre[13];
    ofstream archivo[13];
    ostringstream flujo[13];
    //
    vector <vector<double> > entropia_2D;
    vector <double> entropia_1D;
    vector <int> histograma_1D;
    vector <vector<int> > histograma_2D;
    vector <vector<double> > probabilidad_condicionada;
    vector <vector<int> > contador_de_parametro_de_orden;
    vector <vector<double> > parametro_de_orden_microcanonico;
    vector <vector<double> > parametro_de_orden_microcanonico_promedio;
    vector <int> contador_de_parametro_de_orden_1D;
    vector <double> parametro_de_orden_microcanonico_1D;
    vector <double> parametro_de_orden_microcanonico_promedio_1D;
        
    vector <vector<double> > resultados;
    //ARCHIVOS
    ostringstream graf, anim, cad_datos_del_sistema;
    //
    public:
    Wang_Landau_2D(vector<Modelo*> mod, string descripcion, int Caminantes, int l_de_red, double h_exter, double factor_min, char *lim_de_bloques_de_sweeps, int vecinos, int frec_de_anim=0):h(h_exter) {
        especificacion = descripcion;
        f_minima = factor_min;
        frecuencia_de_animacion = frec_de_anim;
        tiempo_de_simulacion = 0;
        contador_de_sweeps = 0;
        contador_de_bloques = 0;
        cantidad_max_de_bloques = string(lim_de_bloques_de_sweeps);
        if (string(lim_de_bloques_de_sweeps) !="no") {
            num_max_de_bloques = atoi(lim_de_bloques_de_sweeps);
            
        }
        contador_de_animacion = 0;
        caminantes = Caminantes;
        M.resize(caminantes);
        for (int i=0; i<caminantes; i++) {
            M[i] = mod[i];
        }
        N = M[0]->numero_de_espines();
        logaritmo_de_factor_de_modificacion = 1;
        L = l_de_red;
        N = L*L;
        condicion_de_histograma = 10;
        f = 1;
        numero_de_refinamientos = 0;
        temperatura = 0.5;
        cantidad_de_datos = 50;
        J1 = M[0]->valor_de_J1();
        J2 = M[0]->valor_de_J2();
        paso=4;
        v = vecinos/2; //no está en uso
        E_min = (vecinos/2)*N;
        s = E_min/4;
        t = vecinos*N/paso + 1;
        if (frecuencia_de_animacion != 0) {
            anima = 1;
        }
        else {
            anima = 0;
        }
        
                //VECTORES
        entropia_1D.resize(t);
        histograma_1D.resize(t);
        entropia_2D.resize(t);
        histograma_2D.resize(t);
        probabilidad_condicionada.resize(t);
        contador_de_parametro_de_orden.resize(t);
        parametro_de_orden_microcanonico.resize(t);
        parametro_de_orden_microcanonico_promedio.resize(t);
        contador_de_parametro_de_orden_1D.resize(t);
        parametro_de_orden_microcanonico_1D.resize(t);
        parametro_de_orden_microcanonico_promedio_1D.resize(t);
        //NOMBRES
        nombre[0]="Histograma_animado";
        nombre[1]="Histograma_1D";
        nombre[2]="Histograma_2D";
        nombre[3]="Entropia_1D";
        nombre[4]="Entropia_2D";
        nombre[5]="Entropia_1D_final";
        nombre[6]="Entropia_2D_final";
        nombre[7]="Parametro_de_orden";
        nombre[8]="Energia";
        nombre[9]="Magnetizacion";
        nombre[10]="Capacidad_calorifica";
        nombre[11]="Susceptibilidad_magnetica";
        nombre[12]="Cumulante_de_binder";
        //CONSTRUCCIÓN DE VECTORES 2D
        for (int i=0;i<t;i++){
            histograma_2D[i].resize(N+1);
            entropia_2D[i].resize(N+1);
            probabilidad_condicionada[i].resize(N+1); // no está en wl2d
            contador_de_parametro_de_orden[i].resize(N+1);
            parametro_de_orden_microcanonico[i].resize(N+1);
        parametro_de_orden_microcanonico_promedio[i].resize(N+1);
        }
        //
        resultados.resize(5);
        for (int i=0;i<5;i++){
            resultados[i].resize(cantidad_de_datos);
        }
        //INICIALIZANDO VECTORES
        reiniciar_entropia_1D();
        reiniciar_entropia_2D();
        reiniciar_histograma_1D(); //no aparece en WL2D
        reiniciar_histograma_2D();
        reiniciar_parametro_de_orden_microcanonico();
        reiniciar_parametro_de_orden_microcanonico_promedio();
        reiniciar_contador_de_parametro_de_orden();
        reiniciar_probabilidad_condicionada(); // no aparece en WL2D
        //ABRIR PIPE A GNUPLOT
        iniciar_conexion_a_gnuplot();
        //CREAR Y ABRIR ARCHIVOS PARA GUARDAR RESULTADOS
        crear_archivos_de_datos();
        //DATOS INICIALES DE CAMINANTES
        for (int i=0; i<caminantes; i++) {
            M[i]->Energia = M[i]->energia_configuracional();
            M[i]->Magnetizacion = M[i]->magnetizacion();
            M[i]->Parametro_de_orden = M[i]->parametro_de_orden();
        }
        for (int i=0; i<caminantes; i++) {
            cout << "E" << i << ": " <<M[i]->Energia << "\t" << "M" << i << ": " << M[i]->Magnetizacion << "\t" << "PO" << i << ": " << M[i]->Parametro_de_orden << endl;
        }
        //
        estimar_entropia_2D();
        //
        cout << "Tiempo total de simulación: " << tiempo_de_simulacion << endl;
        cout << "Número de bloques totales: " << contador_de_bloques << endl;
        //
        calcular_parametro_de_orden_microcanonico_promedio();
        //
        estimar_entropia_1D();
        //
        entropia_1D_a_archivo();
        //
        entropia_2D_a_archivo();
        //
        //imprimir_resultados();
        //
        graficar();
    }


    void crear_archivos_de_datos(){
        for (int i=0; i<7; i++) {
            flujo[i].str("");
            flujo[i] << nombre[i] << "/" << nombre[i] << "_" << especificacion << ".txt";
            estatus=mkdir(nombre[i].c_str(),S_IRWXU);
            if (!estatus) {
                cout << "Carpeta " << nombre[i] << " creada." << endl;
            }
            else{
                cout << "Carpeta " << nombre[i] << " ya existe." << endl;
            }
            archivo[i].open(flujo[i].str().c_str());
        }
    }
    
    
    void iniciar_conexion_a_gnuplot(){
        anim.str("");
//        anim << "set term wxt 0" << endl;
        anim << "set xrange [" << -E_min << ":" << E_min << "]" << endl;
        anim << "set title 'Histograma animado'" << endl;
        gp[0] = popen("gnuplot -persist", "w");
        //
        graf.str("");
//        graf << "set term wxt 1 persist" << endl;
        gp[1] = popen("gnuplot -persist", "w");
    }
    
    void grafica_evolucion_de_histograma(){
        anim << "plot \"" << flujo[0].str().c_str() << "\"" << endl;
        fprintf(gp[0], "%s\n", anim.str().c_str());
        anim.str("");
        fflush(gp[0]);
    }
    
    
    void datos_de_histograma(){
        archivo[0].open(flujo[0].str().c_str());
        archivo[0] << "#" <<endl;
        for (int i=0;i<t;i++){
            for (int j=0; j<N+1; j++) {
                if (histograma_2D[i][j]!=-1) {
                archivo[0] << 4*i-E_min << "\t" << 2*j-N << "\t" << histograma_2D[i][j] << endl;
                }
            }
            
        }
        archivo[0].close();
    }
    
    int condicion_de_bloques_max(){
        if (cantidad_max_de_bloques == "no") {
            proposicion_de_bloques_max = true;
        }
        else {
            proposicion_de_bloques_max = contador_de_bloques<num_max_de_bloques;
        }

        return proposicion_de_bloques_max;
    }

    
    void reiniciar_entropia_1D(){
        for (int i=0;i<t;i++){
            entropia_1D[i] = 0;
        }
    }  


    void reiniciar_entropia_2D(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                entropia_2D[i][j] = 0;
            }
        }
    }


    void reiniciar_histograma_1D(){
        for (int i=0;i<t;i++){
            histograma_1D[i] = -1;
        }
    }


    void reiniciar_histograma_2D(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                histograma_2D[i][j] = -1;
            }
        }
    }
    
    
    void reiniciar_probabilidad_condicionada(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                probabilidad_condicionada[i][j] = -1;
            }
        }
    }


    void reiniciar_contador_de_parametro_de_orden(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                contador_de_parametro_de_orden[i][j] = -1;
            }
        }
    }


    void reiniciar_parametro_de_orden_microcanonico(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                parametro_de_orden_microcanonico[i][j] = -1;
            }
        }
    }
    
    
    void reiniciar_parametro_de_orden_microcanonico_promedio(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                parametro_de_orden_microcanonico_promedio[i][j] = -1;
            }
        }
    }
    
    
    void reiniciar_contador_de_parametro_de_orden_1D(){
        for (int i=0;i<t;i++){
            contador_de_parametro_de_orden_1D[i] = -1;
        }
    }
    
    
    void reiniciar_parametro_de_orden_microcanonico_1D(){
        for (int i=0;i<t;i++){
            parametro_de_orden_microcanonico_1D[i] = -1;
        }
    }
    
    
    void reiniciar_parametro_de_orden_microcanonico_promedio_1D(){
        for (int i=0;i<t;i++){
            parametro_de_orden_microcanonico_promedio_1D[i] = -1;
        }
    }




    void alterar_logaritmo_de_factor_de_modificacion(){
        f = f/2;
    }
    
    
    void proponer_cambio_configuracional(Modelo *MOD){
        tiempo_de_simulacion++;
        a = MOD->elegir_espin();
        MOD->Energia_nueva = MOD->Energia + MOD->cambio_en_energia_configuracional(a);
        MOD->Magnetizacion_nueva = MOD->Magnetizacion + MOD->cambio_en_magnetizacion(a);
        MOD->Parametro_de_orden_nuevo = MOD->Parametro_de_orden + MOD->cambio_en_parametro_de_orden(a);
    }
    
    
    void implementar_algoritmo(Modelo *MOD){
        double dS = entropia_2D[(MOD->Energia_nueva)/4 + s][(MOD->Magnetizacion_nueva)/2 + N/2] - entropia_2D[(MOD->Energia)/4 + s][(MOD->Magnetizacion)/2 + N/2];
        double p;
        if (dS <= 0) {
            p=1;
            MOD->aceptar_cambio(a);
            MOD->actualizar_variables();
        }
        else {
            p = exp(-dS);
            double r = drand();
            if (r<=p){
                MOD->aceptar_cambio(a);
                MOD->actualizar_variables();
            }
            else{
            }
        }
        int k = (MOD->Energia)/4 + s;
        int l = (MOD->Magnetizacion)/2 + N/2;
        entropia_2D[k][l] += f;
        histograma_2D[k][l]++;
        contador_de_parametro_de_orden[k][l]++;
        parametro_de_orden_microcanonico[k][l] +=  abs(MOD->Parametro_de_orden);
    }
    
    
    void sweep(){
        contador_de_sweeps++;
        for (int i=0;i<N/caminantes;i++){
            for (int j=0; j<caminantes; j++) {
                proponer_cambio_configuracional(M[j]);
                implementar_algoritmo(M[j]);
            }
        }
    }
    
    
    void varios_sweeps(){
        contador_de_bloques++;
        for (int i=0; i<condicion_de_histograma; i++) {
            sweep();
            contador_de_animacion++;
            if (contador_de_sweeps<frecuencia_de_animacion) {
            }
            else {
                contador_de_animacion = 0;
                if (anima==1) {
                    datos_de_histograma();
                    grafica_evolucion_de_histograma();
                }
            }
        }
    }

    
    void renormalizar_entropia_2D(){
        double min=0;
        int k=0;
        int l=0;
        while (entropia_2D [k][l] == 0){
            l++;
            if (l == N) {
                k++;
                l=0;
            }
        }
        min = entropia_2D[k][l];
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if (entropia_2D[i][j]!=0){
                    if (entropia_2D[i][j] <= min){
                        min = entropia_2D[i][j];
                    }
                }
            }
        }
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if (entropia_2D[i][j] !=0 ){
                    entropia_2D[i][j] -= min-1;
                }
            }
        }
    }
    
//Revisar esta renormalización
    void renormalizar_entropia_1D(){
        double norma=0;
        for (int i=0;i<N;i++){
            if (entropia_1D[i] <= entropia_1D[i+1]){
                norma = entropia_1D[i+1];
            }
        }
        for (int i=0;i<t;i++){
            if (entropia_1D[i]!=0){
                entropia_1D[i] -= norma;
            }
        }
    }



    void estimar_entropia_2D(){
        while (f>f_minima && condicion_de_bloques_max()){
            cout << "f: " << f << endl;
            varios_sweeps();
            int k=0;
            if (!condicion_de_bloques_max()) {
                k=t;
            }
            while (k<t){
                int l=0;
                while (l<N+1){
                    if ((histograma_2D[k][l] > condicion_de_histograma-1) || (histograma_2D[k][l] == -1)) {
                            l++;
                    }
                    else {
                        varios_sweeps();
                        if (!condicion_de_bloques_max()) {
                            l=N+1;
                            k=t;
                        }
                    }
                }
                k++;
            }
            //
            renormalizar_entropia_2D();
/*            archivo[4] << endl << endl;
            archivo[2] << endl << endl;
            for (int i=0; i<t; i++){
                for (int j=0; j<N+1; j++){
                    if (entropia_2D[i][j] != 0) {
                        archivo[4] << 4*i-E_min << "\t" << 2*j-N << "\t" << entropia_2D[i][j] << endl;
                        archivo[2] << 4*i-E_min << "\t" << 2*j-N << "\t" << histograma_2D[i][j] << endl;
                    }
                }
            }*/
            alterar_logaritmo_de_factor_de_modificacion();
            numero_de_refinamientos++;
            reiniciar_histograma_2D();
        }
    }


    void estimar_entropia_1D(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N;j++){
                if (entropia_2D[i][j] != 0){
                    entropia_1D[i] += log(entropia_2D[i][j]);
                }
            }
//            archivo[3] << 4*i-E_min << "\t" << entropia_1D[i] << endl;
        }
    }


    void calcular_parametro_de_orden_microcanonico_promedio(){
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if ( contador_de_parametro_de_orden[i][j] != -1){
                    parametro_de_orden_microcanonico_promedio[i][j] = ((parametro_de_orden_microcanonico[i][j] + 1.)/(contador_de_parametro_de_orden[i][j] + 1.));
                }
            }
        }
    }
    
    void entropia_1D_a_archivo(){
        archivo[5] << "Energía \t Magnetización \t Entropía \t P.O." << endl;
        for (int i=0;i<t;i++){
            if (entropia_1D[i]!=-1) {
                archivo[5] << 4*i-E_min << "\t" << entropia_1D[i] << parametro_de_orden_microcanonico_promedio_1D[i] << endl;
            }
        }
        archivo[5].close();
    }
    
    
    void entropia_2D_a_archivo(){
        archivo[6] << "Energía \t Magnetización \t Entropía \t P.O." << endl;
        for (int i=0;i<t;i++){
            for (int j=0;j<N+1;j++){
                if ( entropia_2D[i][j] != 0){
                    //cout << "chechando: " << entropia_2D[i][j] << endl;
                    archivo[6] << 4*i-E_min << "\t" << 2*j-N << "\t" << entropia_2D[i][j] << "\t" << parametro_de_orden_microcanonico_promedio[i][j] << endl;
                }
            }
        }
        archivo[6].close();
    }
    
    
    void encontrar_constante_de_renormalizacion(){
        constante_de_renormalizacion = 0;
        int E, M;
        cout << "valor de beta: " << beta << endl;
        for (int i=0;i<t;i++){
            E = 4*i-E_min;
            for (int j=0;j<N+1;j++){
                if (entropia_2D[i][j]!=0){
                    M = 2*j-N;
                    if ((entropia_2D[i][j] - beta*(E-J2*h*M))>=constante_de_renormalizacion){
                constante_de_renormalizacion = entropia_2D[i][j] - beta*(E-J2*h*M);
                    }
                }
            }
        }
    }
    
    
    void obtener_variables_termodinamicas(int k){
        double funcion_de_particion = 0;
        double funcion_de_particion_2D = 0;
        double energia_promedio = 0;
        double energia_promedio_cuadrada = 0;
        double magnetizacion_promedio = 0;
        double magnetizacion_cuadrada_promedio = 0;
        double magnetizacion_cuarta_promedio = 0;
        double parametro_de_orden_promedio = 0;
        double exponente=0;
        double E=0;
        double M=0;
        for (int i=0;i<t;i++){
            E = 4*i-E_min;
            for (int j=0;j<N+1;j++){
                if(parametro_de_orden_microcanonico_promedio[i][j] != -1){ //aquí no utilizo prob condic
                    M = 2*j-N;
                    exponente = entropia_2D[i][j] - (beta)*(E-J2*h*M) - constante_de_renormalizacion;
                    funcion_de_particion_2D += exp(exponente);
                    energia_promedio += (abs(E))*exp(exponente);
                    energia_promedio_cuadrada += (abs(E*E))*exp(exponente);
                    magnetizacion_promedio += (abs(M))*exp(exponente);
                    magnetizacion_cuadrada_promedio += (abs(M*M))*exp(exponente);
                    magnetizacion_cuarta_promedio += pow(M,4.)*exp(exponente);
                    parametro_de_orden_promedio += (parametro_de_orden_microcanonico_promedio[i][j])*exp(exponente);
                }
            }
        }
        resultados[0][k] = (energia_promedio_cuadrada/funcion_de_particion_2D - (energia_promedio*energia_promedio)/(funcion_de_particion_2D*funcion_de_particion_2D))*(1./N);
        resultados[1][k] = magnetizacion_promedio/(N*funcion_de_particion_2D);
        resultados[2][k] = (magnetizacion_cuadrada_promedio/funcion_de_particion_2D - (magnetizacion_promedio*magnetizacion_promedio)/(funcion_de_particion_2D*funcion_de_particion_2D))/N;
        resultados[3][k] = 1 - ((magnetizacion_cuarta_promedio*funcion_de_particion_2D)/(3*magnetizacion_cuadrada_promedio*magnetizacion_cuadrada_promedio));
        resultados[4][k] = parametro_de_orden_promedio/(N*funcion_de_particion_2D);
    }
    
    
/*    void imprimir_resultados(){
        arch_cap_cal << endl << endl;
        arch_magnet << endl << endl;
        arch_sucep_magnet << endl << endl;
        arch_cbinder << endl << endl;
        arch_po << endl << endl;
        for (int i=0;i<cantidad_de_datos;i++){
            beta = 1/temperatura;
            encontrar_constante_de_renormalizacion();
            obtener_variables_termodinamicas(i);
            arch_cap_cal << temperatura << "\t" << resultados[0][i];
            arch_cap_cal << endl;
//            cout << temperatura << "\t" << resultados[0][i];
//            cout << endl;
            arch_magnet << temperatura << "\t" << resultados[1][i];
            arch_magnet << endl;
            arch_sucep_magnet << temperatura << "\t" << resultados[2][i];
            arch_sucep_magnet << endl;
            arch_cbinder << temperatura << "\t" << resultados[3][i];
            arch_cbinder << endl;
            arch_po << temperatura << "\t" << resultados[4][i];
            arch_po << endl;
            temperatura += 0.1;
        }
    //  cout << "número de refinamientos: " << numero_de_refinamientos -1 << endl;
    }*/
    
    
    void graficar(){
        graf << "encabezado = system('head -1 " << flujo[6].str().c_str() << " ')" << endl;
        graf << "set xlabel word(encabezado,1)" << endl;
        graf << "set ylabel word(encabezado,2)" << endl;
        graf << "set zlabel word(encabezado,3)" << endl;
        graf << "splot \"" << flujo[6].str().c_str() << "\" using 1:2:3 title columnheader" << endl;
        fprintf(gp[1], "%s\n", graf.str().c_str());
        fflush(gp[1]);
/*        cout << endl;
        int espera;
        cin >> espera;*/
    }
    
};


