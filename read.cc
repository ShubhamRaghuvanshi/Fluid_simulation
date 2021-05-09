#include <GL/glut.h>
#include"Fluid_class.h" 
#include"Geometry.h" 
#include <iostream>
#include<fstream>
#include <math.h>

using namespace std;


//1000/60

const float framerate = 1000/60; //milisecond

const int n=256, m=256, scale=3;

int grid_size = n*m;


int i_frame=0, n_frame=1000;

vector<float> Vx, Vy;
vector<int> Bx, By;

float t=0;
float tb=0,te=50, vx, vy, x,y, bx, by, p ;


void init(){
  glClearColor(0,0,0,1);
}




ifstream velocity;
void fillvelocity(){


//  while(!velocity.eof()) {
    
//      velocity>>t>>x>>y>>vx>>vy;
 
      
       
 //     if(t == tp){
     
   //   cout<<"match "<<t<<"  "<<tp<<endl;
 
      
  //      Vx.push_back(vx);
   //     Vy.push_back(vy);
      
        for(int i=0; i<n*m; i++){
          velocity>>i_frame>>x>>y>>vx>>vy>>p>>t;
          Vx.push_back(vx);
          Vy.push_back(vy);
        }
      
    //   break;
    //  }                
//  }
  
}


void display(){ 
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  
  Draw_points( Bx, By);
  
  if(i_frame < n_frame){
    fillvelocity();

    Draw_2dVector_Field( Vx, Vy);

    cout<<i_frame<<" "<<" "<<t<<" "<<Vx.size()<<"  "<<Vy.size()<<endl;


    if(i_frame<  n_frame -1 ){
      Vx.clear();
      Vy.clear();  
    }
    i_frame++;  
  }
  else{

    Draw_2dVector_Field( Vx, Vy);
  }
  

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

  ifstream boundary;

  boundary.open("boundary_pipe.txt");
  velocity.open("fluid_data_pipe.txt");

 
  while(!boundary.eof()) {
  
    boundary>>bx>>by;
    Bx.push_back(bx);
    By.push_back(by);
  }
  
  boundary.close();

  
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































