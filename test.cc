#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>



using namespace std;


int NG=11, MG=10 ,SCALE = 4;
int SG = NG*MG;
float h=1.0/float(NG-1);

int Index(int x, int y){
  return x + NG*y  ;
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


// divergence of vectr field f at point x y 
vector<float> Divergence(vector<float> fx, vector<float> fy, int schm){

  vector<float> dfx, dfy, div;
  
  dfx  = d_dx(fx, schm);
  dfy  = d_dy(fy, schm);

  for(int i=0; i<NG*NG; i++){
    div.push_back(dfx[i] + dfy[i]);
  }
  
  return div;
}


int main(){

  cout<<h<<endl<<endl;

  vector<float> fx, fy, div, dx, dy;
  
  float x,y;
  
  for(int j=0; j<NG; j++){
    for(int i=0; i<NG; i++){
      x = float(i)*h;
      y = float(j)*h; 
      fx.push_back( x );  
      fy.push_back( y );
    }
  }

  div  = Divergence(fx,fy, 0);
 // dx = d_dy(fx,0);
  
/*  
  for(int j=0; j<NG; j++){
    for(int i=0; i<NG; i++){
    
      cout<<setprecision(4)<<setw(12)<<fx[Index(i,j)];
    
    }  cout<<endl;
  }
*/
  
  cout<<endl<<"******************************"<<endl;

  for(int j=0; j<NG; j++){
    for(int i=0; i<NG; i++){
        
      cout<<setprecision(4)<<setw(12)<<div[Index(i,j)];
    }  cout<<endl;
  }
  
  
  return 0;
}




























