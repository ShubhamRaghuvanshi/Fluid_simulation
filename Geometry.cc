#include <GL/glut.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include<fstream>



//#include "Geometry.h"

using namespace std;

//#define float float 

const int NG=256, MG=256 , SCALE = 3;

float qmin = 1e-15;

//Adds points
void point(int x , int y, vector<int> &bx, vector<int> &by ){
  bx.push_back( round(x) );
  by.push_back( round(y) );
}

/*
//plot of fx
void plot(float f[NG]){
  glColor3f(1,1,1);
  glBegin(GL_POINTS);
    for(int i=0; i<NG; i++){
    	if( f[i]>=Y1 && f[i]<=Y2 )
      	glVertex2f(x, f[i]);            
			x = x1 + (i+1)*DX;
			if(abs(x) < DX) glVertex2f(x, Y1);
    }
  glEnd();  
}
*/


//Draws white points
void Draw_points(vector<int> px, vector<int> py){

  glPointSize(SCALE); 
  glColor3f(1,1,1);

  if(px.size() == py.size()){
    glBegin(GL_POINTS);
      for(int i=0; i<px.size(); i++){
        glVertex2f(px[i], py[i]);            
      }
    glEnd();  
  }
  else{
    cout<<"Error in drawing points:: Array sizes not consistent"<<endl;
  }
}

//adds coordinates of pixels along a strainght line between two points
void Line( float x1, float y1, float x2, float y2, vector<int> &bx, vector<int> &by ){

  if( x1>=0 && x1<NG && y1>=0 && y1<NG ){
    float y, m;
    if(x2!=x1){    
      m = (y2-y1)/(x2-x1);     
      for(int x = round(x1); x<=x2; x++){
        bx.push_back(x);
        y =  y1 + m*(x-x1); 
        by.push_back( round(y) );  
      }  
    }
    else {
      for(int y = round(y1); y<=y2; y++){
        bx.push_back(x1);
        by.push_back(y);
      }
    }
  } 
  else{
    cout<<"Error plottig line: End points out of bounds"<<endl;
  }
}


//points on a circle
void Circle( float c_x, float c_y, float r, vector<int> &bx, vector<int> &by){
  float x, y;
  int n = sqrt(2)*M_PI*r;
  float d_theta  = 2.0*M_PI/float(n); 
  for(float theta=0; theta < 2*M_PI; theta = theta + d_theta ){
    x = c_x + r*cos(theta) ;
    y = c_y + r*sin(theta) ;
    bx.push_back(x);
    by.push_back(y);
  }
  
  n = sqrt(2)*M_PI*(r+1);
  d_theta  = 2.0*M_PI/float(n); 
  for(float theta=0; theta < 2*M_PI; theta = theta + d_theta ){
    x = c_x + (r+1)*cos(theta) ;
    y = c_y + (r+1)*sin(theta) ;
    bx.push_back(x);
    by.push_back(y);
  }

  n = sqrt(2)*M_PI*(r-1);
  d_theta  = 2.0*M_PI/float(n); 
  for(float theta=0; theta < 2*M_PI; theta = theta + d_theta ){
    x = c_x + (r-1)*cos(theta) ;
    y = c_y + (r-1)*sin(theta) ;
    bx.push_back(x);
    by.push_back(y);
  }
}

void rect( int x1, int y1, int x2, int y2, vector<int> &bx, vector<int> &by ){
  if(x1>=0 && x1<NG && y1>=0 && y1<NG && x2>=0 && x2<NG && y2>=0 && y2<NG ){
  
    for(int i=x1; i<=x2; i++){
      bx.push_back(i);
      by.push_back(y1);
    }  

    for(int j=y1; j<=y2; j++){
      bx.push_back(x2);
      by.push_back(j);
    }  

    for(int i=x1; i<=x2; i++){
      bx.push_back(i);
      by.push_back(y2);
    }  

    for(int j=y1; j<=y2; j++){
      bx.push_back(x1);
      by.push_back(j);
    }    
  } //
  else{
    cout<<"Error in drawing rectangle: points are out of screen"<<endl;
  }
}

void solid_rect( int x1, int y1, int x2, int y2, vector<int> &bx, vector<int> &by ){
  if(x1>=0 && x1<NG && y1>=0 && y1<NG && x2>=0 && x2<NG && y2>=0 && y2<NG ){
    for(int i=x1; i<=x2; i++){
      for(int j=y1; j<=y2; j++){
        bx.push_back(i);
        by.push_back(j);    
      }    
    }
  }  
  else{
    cout<<"Error in drawing rectangle: points are out of screen"<<endl;
  }
}



//Draws white points
void Draw_density_field(float f[MG][NG], float norm, int x1, int y1, int x2, int y2){

	float c;

	float f_max = 0;	
	for(int j=0; j<MG; j++ ){
		for(int i=0; i<NG; i++ ){
		  if(f_max	<	abs(f[j-y1][i-x1]) )
		    f_max	= abs(f[j-y1][i-x1]);				
		}    
	}


	glPointSize(SCALE);
  glBegin(GL_POINTS);
  for(int j=y1; j<=y2; j++ ){
    for(int i=x1; i<=x2; i++ ){
      if( abs(f[j-y1][i-x1]) > qmin){ 
      	c = f[j-y1][i-x1]/norm;
        if(f[j-y1][i-x1] > 0)
          glColor4f(c, 0, 0, c);
        else 
          glColor4f(0,  0, -c, -c);        
      	glVertex2f(i,j);    
      } //min        
    }  //i   
  } //j
  glEnd();  
}

//Draw scalar field on 2D plane
void Draw_scalar_field(float f[MG][NG], int x1, int y1, int x2, int y2){
  
	float f_max = 0;	
	for(int j=y1; j<=y2; j++ ){
		for(int i=x1; i<=x2; i++ ){
		  if(f_max	<	abs(f[j-y1][i-x1]) )
		    f_max	= abs(f[j-y1][i-x1]);				
		}    
	}
	Draw_density_field(f, f_max/2, x1, y1, x2, y2);
}


//Draws 2D vector field on 2D plane
void Draw_2dVector_Field( float fx[MG][NG], float fy[MG][NG], int nx, int ny, int x1, int y1, int x2, int y2){

  int index, x, y ;
  int gapx = round(	float(NG) / float(nx) );
  int gapy = round(	float(NG) / float(ny)	);
  float temp = 0, f_max, f ;
  
  for(int j=y1; j<=y2;  j++ ){
    for(int i=x1; i<=x2; i++ ){
			if(f_max < fx[j-y1][i-x1]*fx[j-y1][i-x1] + fy[j-y1][i-x1]*fy[j-y1][i-x1] )
				f_max = fx[j-y1][i-x1]*fx[j-y1][i-x1] + fy[j-y1][i-x1]*fy[j-y1][i-x1];
		}
	}	
  f_max = sqrt(f_max);
	if (f_max < qmin) f_max=1;
//	f_max=1;
  glPointSize(5);      
	for(int j=y1; j<=y2 - gapy; j = j + 2*gapy ){
    for(int i=x1; i<=x2 - gapx ; i = i + 2*gapx){
        if( abs(fx[j-y1][i-x1]) > qmin || abs(fy[j-y1][i-x1]) > qmin ){    
          f =  fx[j-y1][i-x1]*fx[j-y1][i-x1] + fy[j-y1][i-x1]*fy[j-y1][i-x1];
          f = sqrt(f);
          x = i+round(gapx*fx[j-y1][i-x1]/f);
          y = j+round(gapy*fy[j-y1][i-x1]/f);
          
          glColor4f(0, 1, 0, 2.0*f/f_max);
          glBegin(GL_LINES);        
          glVertex2f(i,j);        
          glVertex2f(x,y);    
          glEnd();
          
          glBegin(GL_POINTS);
            glVertex2f(x,y);
          glEnd();  
        } // min
      } //i  
    } //j
}

//Integrate f  
float Integrate(float f[MG][NG]){
	float F = 0;
  for(int j=0; j<MG; j++ ){
    for(int i=0; i<NG; i++ ){
      F = F + f[j][i];
		}
	}	
  return F;
}


/*
// divergence of vectr field (fx,fy) 
void Divergence(float fx[MG][NG], float fy[MG][NG], float div[MG][NG], int schm){

  float DfxDx[MG][NG], DfyDy[MG][NG] ;
	if(schm==0){
		ddx_cd(fx, DfxDx);
		ddy_cd(fy, DfyDy);	
		for(int j=0; j<MG; j++){
		  for(int i=0; i<NG; i++){
				div[j][i] = DfxDx[j][i]+DfyDy[j][i];		
			}
		}	
	}
	else if(schm==1){
		ddx_bd(fx, DfxDx);
		ddy_bd(fy, DfyDy);	
		for(int j=0; j<MG; j++){
		  for(int i=0; i<NG; i++){
				div[j][i] = DfxDx[j][i]+DfyDy[j][i];		
			}
		}	
	}
}

// to be continued
void Draw_streamlines( vector<float> fx, vector<float> fy, float px, float py){

  float temp_x =10000 , temp_y =10000 , delt=1;

  int k=0, n_points = 4*NG;

  glColor3f(1,1,1);
  glBegin(GL_LINE_STRIP);
  glVertex2f( round(px) , round(py) );     
  while (k<n_points){


    px = px + fx[  px + NG*py  ]*delt;
    py = py + fy[  px + NG*py ]*delt;
    
    if(temp_x ==  px && temp_y == py) break;
    if(px<0 || py<0 || px>NG-1 || py>MG-1 ) break;
  
    temp_x = px;
    temp_y = py;  

    glVertex2f( round(px) , round(py) );
    k++;
    
  }
   
  glEnd();
}  
*/
	
	 
	 























