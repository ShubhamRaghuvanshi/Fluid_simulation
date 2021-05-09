#include <GL/glut.h>
#include"Fluid_class.h" 
#include"Geometry.h" 
#include <iostream>
#include<iomanip>
#include <math.h>

using namespace std;

//1000/60

const float framerate = 1000/60; //milisecond

const int n=256, m=256, scale=3;

float t=0;
int i_frame=0;
    
Fluid fluid;

void init(){
  glClearColor(0,0,0,1);
}


void display(){ 
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    
    cout<<setw(15)<<"time :"<<setw(15)<<t<<setw(15)<<i_frame<<"\n";

//display functions

    Draw_points( fluid.bx, fluid.by);
    Draw_2dVector_Field( fluid.vx, fluid.vy);

    t = t +0.01;
    i_frame++;

    
//    Draw_streamlines( fluid.vx, fluid.vy , 30 , m/2);
//    Draw_streamlines( fluid.vx, fluid.vy , 30 , m/2 +5);
//    Draw_streamlines( fluid.vx, fluid.vy , 30 , m/2 -5);

    
   // cout<<fluid.vx[ 3 + m/2*n ]<<endl;
    
   Draw_scalar_field( fluid.P );

//Draw_scalar_field( fluid.dens );

  // fluid.Den_evolve();
 

  fluid.Evolve();


    glutSwapBuffers();

}


void reshape(int w, int h){

  glViewport(0,0,(GLsizei) w, (GLsizei) h );
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,n,0,m);
  glMatrixMode(GL_MODELVIEW);

}


void timer(int){
    glutPostRedisplay();

//     fluid.Evolve();

    glutTimerFunc(framerate,timer, 0);
}


int main(int argc, char** argv)
{

 fluid.initialize();
 
  
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   
  glutInitWindowPosition(20, 20);
  glutInitWindowSize(n*scale, m*scale );
  glutCreateWindow("Hello world!");

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable( GL_BLEND );
//   glClearColor(1.0,1.0,0.0,0.0);
  

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutTimerFunc(0,timer,0);    
  
  init();
  glutMainLoop();


    return 0;
}































