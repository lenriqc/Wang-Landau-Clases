#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "red.h"

using namespace std;



class Modelo{

  
  public:
    
    double Energia, Energia_nueva;
    double Magnetizacion, Magnetizacion_nueva;
    double Parametro_de_orden, Parametro_de_orden_nuevo;

    virtual double valor_de_h()=0;
//    virtual void aumentar_h(double dx)=0;
//    virtual void disminuir_h(double dx)=0;

    virtual int energia_configuracional()=0;
    virtual int magnetizacion()=0;
    virtual int parametro_de_orden()=0;
    virtual int cambio_en_parametro_de_orden(int a)=0;
    virtual int cambio_en_energia_configuracional(int a)=0;
    virtual int cambio_en_magnetizacion(int a)=0;
    
    virtual void aceptar_cambio(int a)=0;
    virtual int elegir_espin()=0;
    virtual int numero_de_espines()=0;
    virtual void actualizar_variables()=0;
    
    virtual int valor_de_J1()=0;
    virtual int valor_de_J2()=0;


};



class Ising: public Modelo{
  
  private:
    double h;
    Red* R;
    int N;
    int L;
    int J1,J2;
    string tipo_de_modelo;

  public:
    Ising(Red* r, int H){
        R=r;
        h = H;
        N = R->numero_de_espines();
        L = R->numero_de_espines_por_lado();
        J1 = 1;
        J2 = -1;
        tipo_de_modelo = "Ising-ferro";
    }
    

    int valor_de_J1(){
        return J1;
    }
    
    
    int valor_de_J2(){
        return J2;
    }
        
    
    double valor_de_h(){
        return h;
    }
    
    
    int energia_configuracional(){  
        double E=0;
        for(int i=0; i<N; i++){
            E -= 0.5*(R->valor_del_espin(i))*(R->suma_de_vecinos(i));
        }
        return (J1*E);
    }

      
    int magnetizacion(){
      int M=0;
      for(int i=0; i<N; i++){
          M += R->valor_del_espin(i);
      }
      return (M);
    }

    int parametro_de_orden(){
      int M=0;
      for(int i=0; i<N; i++){
          M += R->valor_del_espin(i);
      }
      return (M);
    }
  
    int cambio_en_parametro_de_orden(int a){
        int dM = -2*(R->valor_del_espin(a));
        return (dM);
    }
 
    
    
    int cambio_en_energia_configuracional(int a){
      int dE = J1*2*(R->valor_del_espin(a))*(R->suma_de_vecinos(a));
      return (dE);
    }
    
    
    int cambio_en_magnetizacion(int a){
      int dM = -2*(R->valor_del_espin(a));
      return (dM);
    }
    
    
    int cambio_en_magnetizacion_staggered(int a){
      int dMS;
      if (R->valor_de_red_reordenada(a) < (N/2)){
	dMS = -2*(R->valor_del_espin(a));
      }
      else{
	dMS = 2*(R->valor_del_espin(a));
      }
      return (dMS);
    }
    
    
    int numero_de_espines(){
        return N;
    }
    
    
    int elegir_espin(){
        return (intrand(R->numero_de_espines()));
    }
    
    
    void aceptar_cambio(int a){
        R->cambiar_valor_de_espin(a);
    }
    
    
    void actualizar_variables(){
        Energia = Energia_nueva;
        Magnetizacion = Magnetizacion_nueva;
        Parametro_de_orden = Parametro_de_orden_nuevo;
    }
    
};





class Antiferro: public Modelo{
    //Falta declarar funciones que ven valores a variables que guarden "energía", "energía_nueva"
  
  private:
    double h;
    Red* R;
    int N;
    int L;
    int J1;
    int J2;
    int forma_de_red;
    string tipo_de_modelo;

  public:
    Antiferro(Red* r, int H,int Tipo_de_red){
        R=r;
        h = H;
        forma_de_red = Tipo_de_red;
        N = R->numero_de_espines();
        L = R->numero_de_espines_por_lado();
        J1=-1;
        J2=-1;
        tipo_de_modelo = "Ising-antiferro";
        Energia = energia_configuracional();
        //Energia_nueva = Energia + cambio_en_energia_configuracional();

    }
    
    int valor_de_J1(){
      return J1;
    }
    
    
    int valor_de_J2(){
      return J2;
    }
    
    
    double valor_de_h(){
      return h;
    }
    

    int energia_configuracional(){  
        double E=0;
            for(int i=0; i<N; i++){
                E -= 0.5*(R->valor_del_espin(i))*(R->suma_de_vecinos(i));
            }
        return (J1*E);
    }


    int magnetizacion(){
        int M=0;
            for(int i=0; i<N; i++){
                M += R->valor_del_espin(i);
            }
        return (M);
    }

 

    int parametro_de_orden(){
    
        int po;
      
        if(forma_de_red == 1){
            int M1=0;
            int M2=0;
            for (int i=0; i<L;i+=2){
                for (int j=0;j<L;j+=2){
                    M1 += R->valor_del_espin(L*i+j);
                }
            }
            for (int i=1; i<L;i+=2){
                for (int j=1;j<L;j+=2){
                    M1 += R->valor_del_espin(L*i+j);
                }
            }
            for (int i=0; i<L;i+=2){
                for (int j=1;j<L;j+=2){
                    M2 += R->valor_del_espin(L*i+j);
                }
            }
            for (int i=1; i<L;i+=2){
                for (int j=0;j<L;j+=2){
                    M2 += R->valor_del_espin(L*i+j);
                }
            }
            po = M1-M2;
        }
        
        if(forma_de_red == 2){
            int M1=0;
            int M2=0;
            int M3=0;
            int M=0;
            for(int i=0; i<L; i+=3){
                for(int j=0; j<L; j+=3){
                    M1 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=1; i<L; i+=3){
                for(int j=2; j<L; j+=3){
                    M1 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=2; i<L; i+=3){
                for(int j=1; j<L; j+=3){
                    M1 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=0; i<L; i+=3){
                for(int j=1; j<L; j+=3){
                    M2 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=1; i<L; i+=3){
                for(int j=0; j<L; j+=3){
                    M2 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=2; i<L; i+=3){
                for(int j=2; j<L; j+=3){
                    M2 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=0; i<L; i+=3){
                for(int j=2; j<L; j+=3){
                    M3 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=1; i<L; i+=3){
                for(int j=1; j<L; j+=3){
                    M3 += R->valor_del_espin(j+L*i);
                }
            }
            for(int i=2; i<L; i+=3){
                for(int j=0; j<L; j+=3){
                    M3 += R->valor_del_espin(j+L*i);
                }
            }
            M1 = abs(M1);
            M2 = abs(M2);
            M3 = abs(M3);
            po = M1+M2+M3;
        }
        return(po);
    }
    
    
    int cambio_en_parametro_de_orden(int a){
        int dMS;
        if(forma_de_red == 1){
            if (R->valor_de_red_reordenada(a) < (N/2)){
                dMS = -2*(R->valor_del_espin(a));
            }
            else{
                dMS = 2*(R->valor_del_espin(a));
            }
        }
        if(forma_de_red == 2){ //ojo aquí tal vez sólo necesito un if y un else
            if (R->valor_de_red_reordenada(a) < (N/3)){
                dMS = -2*(R->valor_del_espin(a));
            }
            if (R->valor_de_red_reordenada(a) > (2*N/3)){
                dMS = 2*(R->valor_del_espin(a));
            }
            else{
                dMS = -2*(R->valor_del_espin(a));
            }
        }
        return (dMS);
    }

      
    
    int cambio_en_energia_configuracional(int a){
      int dE = J1*2*(R->valor_del_espin(a))*(R->suma_de_vecinos(a));
      return (dE);
    }
    
      
    int cambio_en_magnetizacion(int a){
      int dM = -2*(R->valor_del_espin(a));
      return (dM);
    }
    
    
    int cambio_en_magnetizacion_staggered(int a){
      int dMS;      
      if (R->valor_de_red_reordenada(a) < (N/2)){
	dMS = -2*(R->valor_del_espin(a));
      }
      else{
	dMS = 2*(R->valor_del_espin(a));
      }
      return (dMS);
    }
    
    
    int numero_de_espines()
    {
      return N;
    }
    
    int elegir_espin(){
      return (intrand(R->numero_de_espines()));
    }
    
    void aceptar_cambio(int a){
      R->cambiar_valor_de_espin(a);
    }
    
    void actualizar_variables(){
        Energia = Energia_nueva;
        Magnetizacion = Magnetizacion_nueva;
        Parametro_de_orden = Parametro_de_orden_nuevo;
    }
      
};

