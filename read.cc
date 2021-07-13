#include"../Fluid_class.cc" 
#include "../Geometry.cc"

const float framerate = 1000/30; //milisecond

vector<int> Bx, By;
int Nx=70, Ny=70;
float Vx[MG][NG]= {}, Vy[MG][NG]= {}, p[MG][NG]= {}, w[MG][NG]= {}, psi[MG][NG]= {}; 
  
const int len=100;  
char strng[len] ;        
float t;    

int x1 = NG/2-NG/4, wy1 = MG/2-MG/4;
int x2 = NG/2+NG/4, wy2 = MG/2-MG/4;	
int x3 = NG/2+NG/4, wy3 = MG/2+MG/4;
int x4 = NG/2-NG/4, wy4 = MG/2+MG/4;


void init(){
  glClearColor(0,0,0,1);
}

void display(){ 
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

	glColor3f(1, 1, 1);  
	glRasterPos2i(10, MG);
//	sprintf(strng, "t=%0.3f     Re=%0.3f      Eu=%0.3f      T=%0.3f", t, RE, EU, T);
//  for (int i = 0; i < len; i++) {
 //   glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (int)strng[i]);
//	}


//display functions
//	Draw_points(Bx,	By);
	Draw_2dVector_Field( Vx, Vy, Nx, Ny,  x1, wy1, x3-20, wy3-20);
//	Draw_scalar_field(psi, x1+1, wy1+1, x3-1, wy3-2);
//Draw_scalar_field(w, x1, wy1, x3, wy3);
//	Draw_density_field(d,1);
 // Advect_C(d, Vx, Vy, DT);
 // ConvectDiffuse(d, Vx, Vy, 5, DT);
//	cout<<Integrate(d)<<endl;
//  Draw_scalar_field( fluid.P );
//	T = T0 + it*DT;
//	it++;
	glutSwapBuffers();
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
	char filename[25];
	float EU=1, RE=400;
	sprintf(filename,"flow_%0.3f_%0.3f", EU, RE);
	
	read_from_file(t, Vx, Vy, p, w, psi, Bx, By, filename);
//	for(int j=MG/2-MG/4; j<=MG/2+MG/4; j++){
//		for(int i=20; i<25; i++){
//			d[j][i]=3;
//		}
//	}


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































