// OpenGL y GLUT basicos
// Compilar con la opcion -lglut  para ligar con la libreria de GLUT

//freeglut-dev: hay que instalar esto en UBUNTU

// g++ ddd.cpp -o ddd lglut
 // uso del teclado
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include "algoritmo3.h"

Red* r;
Modelo* m;
Algoritmo* a;


 
using namespace std;

// int N=100;
int L;
double delta_x;


// double puntero_a_red(Red*r){
//   return r
  
// }

void inic_glut(int argc, char** argv, Red*r) {
// 	int N = r->numero_de_espines();

	L=r->numero_de_espines_por_lado();
	  
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(600,600);
	glutInitWindowPosition(500,500);
	glutCreateWindow("ISING");
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-L/2,L,-L/2,L/2,-L/2,L/2);
	glMatrixMode(GL_MODELVIEW);
}



void display() {

  
    glClear(GL_COLOR_BUFFER_BIT);        // limpiar pantalla

	for(int j=0; j<L; j++)
	{
		for(int i=0; i<L; i++)
		{
			glBegin(GL_POLYGON);
			if(r->valor_del_espin(j+L*i)==1){
				glColor3f(1.,0.,0.);
			}
			else{
				glColor3f(0.,0.,1.);
			}
			glVertex2f(-(L - (L-j))/2. + i ,(L/2 - 1) - j);
			glVertex2f(-((L - (L-j))/2. - 1) + i,(L/2 - 1) - j);
			glVertex2f(-((L - (L-j))/2. - 1) + i,L/2 - j);
			glVertex2f(-(L - (L-j))/2. + i,L/2 - j);
	
			glEnd();  
		}
	}

	
	
                                                 // terminar de dibujar el triángulo
    glFlush();                                              // mandar los comandos ya
    
}
 
void actualizar(){
	a->sweep();
	glutPostRedisplay();
	glutSwapBuffers();
}

void teclado(unsigned char tecla, int x, int y){
  switch(tecla){
    case 'q':
      exit(0);
    case 'u':
      a->disminuir_beta(delta_x);
      cout << "T = " << 1/(a->valor_de_beta())<< endl;
      break;
    case 'd':
      a->aumentar_beta(delta_x);
      cout << "T = " << 1/(a->valor_de_beta())<< endl;
      break;
    case '-':
      m->disminuir_h(delta_x);
      cout << "h = " << m->valor_de_h() << endl;
      break;
    case '+':
      m->aumentar_h(delta_x);
      cout << "h = " << m->valor_de_h() << endl;
      break;
  }

}

void raton(int boton, int estado, int x, int y){
  switch (boton){
    case GLUT_LEFT_BUTTON:
      if(estado ==GLUT_DOWN){
	      glutIdleFunc(actualizar);
      }
      break;
      
    case GLUT_RIGHT_BUTTON:
      if(estado ==GLUT_UP){
	      glutIdleFunc(NULL);
      }
      break;
    }
}
 

 
 
int main(int argc, char** argv) {

  srand(time(0));
	
  


    
    int tipo_de_red;
    cout << "Dame el número para escoger el tipo de red: " << endl;
    cout << "1: cuadrada\n2: triangular\n--> ";
    cin >> tipo_de_red;
    
    int tipo_de_modelo;
    cout << "Elige un número para el modelo a simular: " << endl;
    cout << "1: Ising\n2: Antiferromagnético\n--> ";
    cin >> tipo_de_modelo;
    
  
  
//   if (tipo=="cuadrada"){
    cout << "Dame el número de espines por lado: ";
    int LL;
    cin >> LL;

//     cout << LL << endl;
    
    cout << "Dame el valor del campo externo: ";
    double h;
    cin >> h;
    
    cout << "Dame el valor inicial de beta: ";
    double beta;
    cin >> beta;
  
    cout << "Dame el intervalo para cambiar beta: ";
    
    cin >> delta_x;
  
    
    
  switch (tipo_de_red){
  
  case 1:
    r = new Red_Cuadrada(LL);
    break;
    
  case 2:
    r = new Red_Triangular(LL);
    break;
  }
  
    
  
  switch (tipo_de_modelo){
   
    case 1:
      m = new Ising(r,h);
      break;
      
    case 2:
      m = new Antiferro(r,h);
      break;
  }
    
    
  a = new Metropolis(m, beta);
//   }
  
  
  
  
  
  
  
  inic_glut(argc,argv,r);

  glutDisplayFunc(display);  // decirle a glut cual función usar para desplegar 
    
  glutIdleFunc(actualizar);
	
  glutKeyboardFunc(teclado);
	
  glutMouseFunc(raton);
	
  glutMainLoop();                         // entrar en el bucle principal de glut
    
}
 
