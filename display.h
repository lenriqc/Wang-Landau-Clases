#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <gl/glut.h>
#include "algoritmo3.h"



class Configuracion{
  public:
    virtual void Inicializar_Glut(int argc, char** argv)=0;
    virtual void Display()=0;
    
};




class Configuracion_Cuadrada: public Configuracion{
  
  private:
    Red* R;
    int L;
   
  public:
    
    Configuracion_Cuadrada(Red* r){
      R=r;
      L = R->numero_de_espines_por_lado();
    }    
    
    void Inicializar_Glut(int argc, char** argv){

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(600,600);
	glutInitWindowPosition(500,500);
	glutCreateWindow("SIMULACIÓN");
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-L/2.,L/2.,-L/2,L/2,-L/2,L/2);
	glMatrixMode(GL_MODELVIEW);
    }
    
    
    void Display(){
      
      glClear(GL_COLOR_BUFFER_BIT);        // limpiar pantalla

	for(int j=0; j<L; j++)
	{
		for(int i=0; i<L; i++)
		{
			glBegin(GL_POLYGON);
			if(R->valor_del_espin(j+L*i)==1){
				glColor3f(1.,0.,0.);
			}
			else{
				glColor3f(0.,0.,1.);
			}
			glVertex2f(-L/2 + i,(L/2 - 1) - j);                                // especificar un vértice en 2D, que consiste en números flotatntes
			glVertex2f(-(L/2 - 1) + i,(L/2 - 1) - j);
			glVertex2f(-(L/2 - 1) + i,L/2 - j);
			glVertex2f(-L/2 + i,L/2 - j);
	
			glEnd();  
		}
	}

      glFlush();
    
    }
};
    
    
    
class Configuracion_Triangular: public Configuracion{
  
  private:
    Red* R;
    int L;
    double radio;
    double x;
    double longitud;
    double altura;
  public:
    
    Configuracion_Triangular(Red* r){
      R=r;
      L = R->numero_de_espines_por_lado();
      radio = sqrt(3);
      x = radio*sqrt(3);
      
    }
        
    void Inicializar_Glut(int argc, char** argv) {
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(1000,1000);
	glutInitWindowPosition(100,100);
	glutCreateWindow("SIMULACIÓN");
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
// 	glOrtho(-L*x/2,L*x,-L*3*radio/2,L*radio/2, -L/2, L/2);
	
	for (int i=0; i<L; i+=2){
	    altura += radio*2;
	}
	for (int i=1; i<L; i+=2){
	  altura += radio;
	}
	
	longitud = x*L + x*(L-1)/2;

// 	int longitud = (L/2+1)*2*radio+(L/2)*radio;
	glOrtho(-x*(L-1)/2,x*L,-longitud/2,longitud/2,-longitud/2,longitud/2);
	glMatrixMode(GL_MODELVIEW);
    }


    void Display() {  
      
      glClear(GL_COLOR_BUFFER_BIT);// limpiar pantalla
      

	  for(int j=0; j<L; j++)
	  {
		  for(int i=0; i<L; i++)
		  {
			  glBegin(GL_POLYGON);
			  if(R->valor_del_espin(j+L*i)==1){
				  glColor3f(1.,0.,0.);
			  }
			  else{
				  glColor3f(0.,0.,1.);
			  }
			  glVertex2f(-((L - (L-j*x))/2.) + i*x , (altura/2. - radio/2.) - j*3*radio/2.);
			  glVertex2f(-((L - (L-j*x))/2.) + i*x , (altura/2. - 3*radio/2.) - j*3*radio/2.);
			  glVertex2f(-((L - (L-j*x))/2.) + x/2. + i*x , (altura/2. - 2*radio) - j*3*radio/2.);
			  glVertex2f(-((L - (L-j*x))/2.) + x + i*x , (altura/2. - 3*radio/2.) - j*3*radio/2.);
			  glVertex2f(-((L - (L-j*x))/2.) + x + i*x , (altura/2. - radio/2.) - j*3*radio/2.);
			  glVertex2f(-((L - (L-j*x))/2.) + x/2. + i*x  , (altura/2.) - j*3*radio/2.);
	  
			  glEnd();  
		  }
	  }

      glFlush();
      
    }
};