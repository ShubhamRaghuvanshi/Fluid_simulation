#include <GL/glut.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "Geometry.h"

using namespace std;


int NG=256, MG=256 , SCALE = 3;
int SG = NG*MG;
double h=1.0;

int Index(int x, int y){
  return x + NG*y  ;
}


//Draw scalar field on 2D plane
void Draw_scalar_field(vector<float> f){

  
  if(f.size() == SG){

    float f_max = f[0];
    
    for(int i=1;i<SG;i++){
      if(f_max < abs(f[i]) )
        f_max = abs(f[i]);
    }

    for(int i=0;i<SG;i++){      
      f[i] = f[i]/f_max;
    }

    glBegin(GL_POINTS);
    for(int j=0; j<MG; j++ ){
      for(int i=0; i<NG; i++ ){
      
        if(f[Index(i,j)] >= 0)
          glColor4f(1, 0, 0, f[Index(i,j)]);
        else 
          glColor4f(0, 0, 1, -f[Index(i,j)]);
        
        glVertex2f(i,j);
      }     
    }
    glEnd();  
  }  //size if 
  else{
    cout<<"Draw scalar field error: Field not defined over whoe domain"<<f.size()<<endl;
  }
}


//mass of the fluid if f is density
float Mass(vector<float> f){

  if(f.size() == SG){

    float mass = 0;
    
    for(int i=0;i<SG;i++){
        mass = mass + f[i];
    }
    return mass;
  }  //size if 
  else{
    cout<<"Field not defined over whoe domain"<<f.size()<<endl;
    return 0;
  }
}
  





//Draws 2D vector field on 2D plane
void Draw_2dVector_Field( vector<float> fx, vector<float> fy){

  int index, x, y ;

  int gap = round(  float(NG) /70.0 );
  
  float temp = fx[0]*fx[0] + fy[0]*fy[0], f_max, f ;
  
  
  glPointSize(5);   
  if(fx.size() == SG && fy.size() == SG ){


    for(int k=1; k<SG; k++){
      temp = fx[k]*fx[k] + fy[k]*fy[k];
      if (f_max < temp)  f_max = temp;   
 
    }  
    f_max = sqrt(f_max);
 
 //   cout<<"max : "<<f_max<<endl;

 
 
    for(int j=0; j<MG - gap; j = j + 2*gap ){
      for(int i=0; i<NG - gap ; i = i + 2*gap ){
        index = Index(i,j);
        f = fx[index]*fx[index] + fy[index]*fy[index];
        f = sqrt(f);
        x = i+round(gap*fx[index]/f);
        y = j+round(gap*fy[index]/f);
        
       
//        glColor3f(0,  f    , 0);        
        glColor4f(0, 1, 0, f/f_max);

        
        glBegin(GL_LINES);        
        glVertex2f(i,j);        
        glVertex2f(x,y);    
        glEnd();
        
        glBegin(GL_POINTS);
          glVertex2f(x,y);
        glEnd();  

         
      } //i loop 
    } //j loop
    

  } //if
  else{
    cout<<"Error in plotting field: Vector field is not defined over whole domain."<<  fx.size() <<"  "<<fy.size() <<endl;
  } 
}

//d/dx field with periodic boundary 
vector<float> d_dx(vector<float> f,  int schm ){
  
  float der;
  vector<float> dfx;  

  if(schm ==0 ) {

    for(int y=0; y<MG; y++){
      for(int x=0; x<NG; x++){

        if(x>0 && x< NG-1 ){
          der = f[Index(x+1,y)] - f[Index(x-1,y)]; 
        }

        else if(x==0){
          der = f[Index(1,y)] - f[Index(NG-1,y)];          
        }
        else if(x==NG-1){
          der = f[Index(0,y)] - f[Index(NG-2,y)];          
        }
        else{}

        der = der/(2.0*h);                   
        dfx.push_back(der);  
      }
    }  
  } //cd

  else if(schm ==1 ) {

    for(int y=0; y<MG; y++){
      for(int x=0; x<NG; x++){

        if(x>0 ){
          der = f[Index(x,y)] - f[Index(x-1,y)];          
        }

        else if(x==0){
          der = f[Index(0,y)] - f[Index(NG-1,y)];
        } 
        else {}
        
        der = der/h;
        dfx.push_back(der);  
      }
    }  
  } //bd


  else {}

  return dfx;
}

//d/dy field with periodic boundary 
vector<float> d_dy(vector<float> f,  int schm){

  float der;
  vector<float> dfy;  

  if(schm ==0 ) {

    for(int y=0; y<MG; y++){
      for(int x=0; x<NG; x++){

        if(y>0 && y< MG-1 ){
          der = f[Index(x,y+1)] - f[Index(x,y-1)];
                  
        }

        else if(y==0){
          der = f[Index(x,1)] - f[Index(x,MG-1)];          
        }
        
        else if(y==MG-1){
          der = f[Index(x,0)] - f[Index(x,MG-2)];          
        }
        else{}

        der = der/(2.0*h); 
        dfy.push_back(der);  
      }
    }  
  } //cd

  else if(schm ==1){
    for(int y=0; y<MG; y++){
      for(int x=0; x<NG; x++){

        if(y>0){
          der = f[Index(x,y)] - f[Index(x,y-1)];  
        }
        else if (y == 0){
          der = f[Index(x,0)] - f[Index(x,MG-1)];        
        }
        else {}
        
        der = der/h;
        dfy.push_back(der);  
      }
    }  
  }   //bd 
  
  else {}
  
  return dfy;

}


//laplacian field with periodic boundary 
vector<float> laplacian(vector<float> f){
  float ddx,ddy;
  vector<float> ddf;  

  double hh = h*h;

  for(int y=0; y<MG; y++){
    for(int x=0; x<NG; x++){

      if(x>0 && x<NG-1){
        ddx = f[Index(x+1,y)] -2.0*f[Index(x,y)] + f[Index(x-1,y)];  
      }
      else if (x == 0){
        ddx = f[Index(1,y)] -2.0*f[Index(0,y)] + f[Index(NG-1,y)];  
      }
      else if (x == NG-1){
        ddx = f[Index(0,y)] -2.0*f[Index(x,y)] + f[Index(x-1,y)];        
      }   
      else {}

      if(y>0 && y<NG-1){
        ddy = f[Index(x,y+1)] -2.0*f[Index(x,y)] + f[Index(x,y-1)];  
      }
      else if (y == 0){
        ddy = f[Index(x,1)] -2.0*f[Index(x,0)] + f[Index(x,NG-1)];  
      }
      else if (y == NG-1){
        ddy = f[Index(x,0)] -2.0*f[Index(x,y)] + f[Index(x,y-1)];        
      }   
      else {}
 
 
      ddf.push_back( (ddx+ddy)/(hh) );
    }
  }  
  return ddf;

}


// divergence of vectr field f at point x y 
vector<float> Divergence(vector<float> fx, vector<float> fy, int schm){

  vector<float> dfx, dfy, div;
  
  dfx  = d_dx(fx, schm);
  dfy  = d_dy(fy, schm);

  for(int i=0; i<SG; i++){
    div.push_back(dfx[i] + dfy[i]);
  }
  
  return div;
}



void Draw_streamlines( vector<float> fx, vector<float> fy, float px, float py){

  //cout<<px<<" "<<py<<"    "<<fx[px + NG*py]<<"  "<<fy[px + NG*py]<<endl;


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



//Draws white points
void Draw_dens_points(vector<int> px, vector<int> py, vector<float> dens ){

  int index;

  glPointSize(SCALE); 
  if(px.size() == py.size() == dens.size() ){

    glBegin(GL_POINTS);
    for(int i=0; i<px.size(); i++){

      if( dens[i] >1  )  dens[i] ==1;
      if( dens[i] < 0 )  dens[i] ==0; 

      glColor4f(1, 1, 1, dens[i] );                
      glVertex2f(px[i], py[i]); 
       
    }//i
    glEnd();  
  }
  else{
    cout<<"Error in drawing density points:: Array sizes not consistent"<<endl;
  }
}




//Adds points
void point(int x , int y, vector<int> &bx, vector<int> &by ){
  bx.push_back( round(x) );
  by.push_back( round(y) );
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
void Circle( float c_x, float c_y, float r, vector<int> &bx, vector<int> &by, int n ){

  float x, y;
  
  float d_theta  = 2.0*M_PI/float(n);
  

  for(float theta=0; theta < 2*M_PI; theta = theta + d_theta ){
    x = c_x + r*cos(theta) ;
    y = c_y + r*sin(theta) ;

    bx.push_back(round(x));
    by.push_back(round(y));


 
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











































