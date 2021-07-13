#include"../Fluid_class.cc" 
#include "../Geometry.cc"
#include <unistd.h>

using namespace std;

const float framerate = 1000/30	; //milisecond

vector<int> Bx, By;
int Nx=70, Ny=70;
float Vx[MG][NG]= {}, Vy[MG][NG]= {}, w[MG][NG]= {}, p[MG][NG]= {}; 
vector<float> VBX, VBY;


float u=1; 
const int len=50;  
char strng[len] ;        

int it=0;   
Fluid fluid;

float t0 = fluid.t;
void init(){
  glClearColor(0,0,0,1);
}


void display(){ 
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

	glColor3f(1, 1, 1);  
	glRasterPos2i(10, MG);
	sprintf(strng, "t=%f     Re=%f      Eu=%f",fluid.t,fluid.Re,fluid.Eu);
  for (int i = 0; i < len; i++) {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (int)strng[i]);
	}

//display functions
	Draw_points(fluid.bx,	fluid.by);

	Draw_2dVector_Field( fluid.vx, fluid.vy, Nx, Ny, 0, 0, NG, MG);
	Draw_scalar_field(fluid.w,0, 0, NG, MG);
//	Draw_scalar_field(f.luid.psi, x1+1, wy1+1, x3-1, wy3-2);
//	Draw_density_field(fluid.psi, 1, x1+1, wy1+1, x3-1, wy3-2);

 //
	//sleep(2);
// to evolve 
//  fluid.Evolve();

	fluid.vorticity_stream();
//	fluid.pressure_corrector();
  fluid.t =t0 + (it+1)*fluid.dt;
	it++;
	glutSwapBuffers();	

/*
	if( it < 500 ){

	glReadPixels(0, 0, NG, MG, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
	fwrite(buffer, 3*NG*MG, 1, ffmpeg);
	}
	else if (it==500){
		pclose(ffmpeg);
		 glutLeaveMainLoop();
	}
	else{}
*/
}


void reshape(int w, int h){
  glViewport(0,0,(GLsizei) w, (GLsizei) h );
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,NG,0,MG+5);
  glMatrixMode(GL_MODELVIEW);
}


void timer(int){
    glutPostRedisplay();
    glutTimerFunc(framerate,timer, 0);
}


int main(int argc, char** argv)
{

	for(int j=0; j<M; j++ ){
		for(int i=0; i<N; i++ ){
			Vx[j][i] = ((float) rand() / (RAND_MAX))/2 ;
			Vy[j][i] =  ((float) rand() / (RAND_MAX))/2 ;

		}	
	}
	
	fluid.initialize(Vx,Vy,p,Bx,By,VBX,VBY);



  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   
  glutInitWindowPosition(20, 20);
  glutInitWindowSize(NG*SCALE, MG*SCALE );
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































