#include <GL/glut.h>
#include"Fluid_class.h" 
#include"Geometry.h" 
#include <iostream>
#include<iomanip>
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

float t=0, dt=0.05;
float tb=0,te=50, vx, vy, x,y, bx, by, p ;


void init(){
  glClearColor(0,0,0,1);
}


vector<float> dens,temp_dens;
void init_density(){
 int index;
 
 for(int i=0; i<grid_size; i++){
  dens.push_back(0);
 }
 
 for( int y=m/2-m/4 ; y<=m/2+m/4   ;y = y+1 ){
//   index = IX( 1,y);
   dens[10 + n*y ] =10;
  }      
}



void Den_evolve(){

/*
  float mass = 0;  
  for(int i=0;i<grid_size;i++){
      mass = mass + dens[i];
  }
  cout<<"Mass = "<<mass<<endl; 
*/
  
  vector<float> rhovx, rhovy, div_rhov; 

  for(int i=0;i<grid_size;i++){
    rhovx.push_back(dens[i]* Vx[i] );
    rhovy.push_back(dens[i]* Vy[i] );    
  }
  
  div_rhov = d_dx(rhovx,1);
  
  float dv=0;
  for(int i=0;i<grid_size;i++){
    dv = dv + abs(div_rhov[i]);
  }  
//  cout<<"div : "<<dv<<endl;
  

  temp_dens = dens;  
  for(int i=0;i<grid_size;i++){
    temp_dens[i] = dens[i] - dt*div_rhov[i];
  }
  
  dens = temp_dens;
}



void display(){ 
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

  Draw_points( Bx, By);
  Draw_2dVector_Field( Vx, Vy);
  

  Draw_scalar_field(dens);

 if(i_frame < 750){
  Den_evolve();  
  tb = tb+dt;
  cout<<setw(15)<<"t :"<<setw(15)<<tb<<setw(15)<<i_frame<<endl;  

  }
  i_frame++;


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

  ifstream boundary, velocity;

  boundary.open("boundary_pipe.txt");
  velocity.open("instantanous_fluid_data_pipe.txt");


  while(!boundary.eof()) {
  
    boundary>>bx>>by;
    Bx.push_back(bx);
    By.push_back(by);
  }
  
       
  for(int i=0; i<n*m; i++){
    velocity>>n_frame>>x>>y>>vx>>vy>>p>>t;
    Vx.push_back(vx);
    Vy.push_back(vy);
  }
         
  init_density();        
         
  boundary.close();
  velocity.close();

  
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































