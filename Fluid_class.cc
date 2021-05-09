
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include<fstream>
#include "Fluid_class.h"
#include "Geometry.h"

using namespace std;



int Fluid::IX(int x, int y){
  return x + N*y  ;
}


Fluid::Fluid(){

   N=256;
   M=256;
   grid_size = N*M;


  for(int i=0; i<grid_size; i++ ){

    temp_dens.push_back(0);
    temp_vx.push_back(0);
    temp_vy.push_back(0); 

    dens.push_back(0);
    vx.push_back(0);
    vy.push_back(0);   
    
//    temp_w.push_back(0);   
    P.push_back(0);   

  }
//  cout<<"init : "<<dens.size()<<" "<<grid_size<<endl;
} //constructor




Fluid::~Fluid(){}; //distructor


void Fluid::init_boundary(){

   //pipe
//   Line(0, M/2+M/4,   N-1, M/2+M/4,   bx, by);
//   Line(0, M/2+M/4+1, N-1, M/2+M/4+1, bx, by); 
//   Line(0, M/2-M/4,   N-1, M/2-M/4,   bx, by);
//   Line(0, M/2-M/4-1, N-1, M/2-M/4-1, bx, by); 


   //boundary
//   Line(1,0,1,M-1, bx, by);
//   Line(0,0,0,M-1, bx, by); 
//   Line(N-1,0,N-1,M-1, bx, by);
//   Line(N-2,0,N-2,M-1, bx, by); 

   Line(0, 7*M/8+1 ,   N-1, 7*M/8+1,   bx, by);
   Line(0, 7*M/8, N-1, 7*M/8, bx, by); 
   Line(0, M/8,   N-1, M/8,   bx, by);
   Line(0, M/8-1, N-1, M/8-1, bx, by); 



  solid_rect( 50, N/2-20, 90, N/2+20, bx, by );


   
   
}


void Fluid::init_density(){

 int index;
/* 
 for(int y=0; y<M; y++){
   for(int x=0; x<N; x++){
    index = IX( x,y);
   
    if(x <= N/4) dens[index] = 1;
    else if( x> N/4 && x <= N/2) dens[index] = 0.5;
    else if( x> N/2 && x <= 3*N/4) dens[index] = -0.5;
    else if( x> 3*N/4 && x <= N-1) dens[index] = -1;
    else{}
   } 
 }
*/           
 for( int y=M/2-M/4 ; y<=M/2+M/4   ;y = y+1 ){
   index = IX( 1,y);
   dens[index] =10;
  }      
}

/*
  for(int y=0 ; y<M; y++){
   for(int x=0; x<N; x++){
      index = IX(x,y);      
      
    if(x <= N/4) { vx[index] = 1; vy[index] = 0;  }
    else if( x> N/4 && x <= N/2) { vx[index] = 1; vy[index] = 1;  }
    else if( x> N/2 && x <= 3*N/4) { vx[index] = 1; vy[index] = -1;  }
    else if( x> 3*N/4 && x <= N-1) { vx[index] = 0; vy[index] = -1;  }
    else{}

*/


void Fluid::init_velocity(){

  int index;
  
  float u=5;

  for(int y=M/8+1   ; y<=7*M/8-1 ; y++){
   for(int x=0; x<N-1 ; x++){  
      index = IX(x,y);
      vx[index]= u;
      vy[index]= 0;  
      }    
  }  
 
  Re = rho*u*del_x/eta;  
  
  if( dt <= 2.0*eta/(u*u*rho) && del_x >= u*dt ){
    cout<<"Stable solution expected  "<<endl;
  }
  else{
    cout<<"Error : Solution might be unstable"<<endl;
  }
  cout<<"Allowed values of grid and time spacing for stable solution are:  dt < "<<2.0*eta/(u*u*rho)<<"  dx >"<<2.0*eta/(u*rho)<<endl<<endl;
  
  cout<<"Reynolds number of flow per unit cell is "<<Re<<endl;

/*
  Re = rho*del_x/eta;
  float RE=0;
  for(int i=0; i<grid_size; i++){
    RE = RE + Re*vx[i]  ;
  }
  cout<<"Overall Reynolds number of the flow is : "<<RE<<endl;
*/

  cout<<"Domain of the analysis is "<<0<<" <= x <= "<<(N-1)*del_x;
  cout<<"  , "<<0<<" <= y <= "<<(M-1)*del_x<<endl;


}

/*
void Fluid::set_vorticity(){

  int index;
  
  for(int y=0; y<N; y++){
    for(int x=0; x<N; x++){
      index = IX(x,y);      
      w[index]=0;             
   }    
  }  
}
*/

////
void Fluid::noslip(){

  for(int i=0; i<by.size(); i++ ){
    vx[IX(bx[i],by[i])] =0;
    vy[IX(bx[i],by[i])] =0;   
  }  

}




void Fluid::initialize(){

  init_boundary();  
 // init_density();
  init_velocity();  
  

  noslip();
}


//void Fluid::solve_steady_incompressible(){

 // noslip();
 
 
//}





vector<float> Fluid::Advect(vector<float> f ){

 float adx, ady ;
 int index;
 vector<float> ad; 

  for(int y=0; y<M; y++){
    for(int x=0; x<N; x++){

      index = IX(x,y);

      if(x>0){
        adx = vx[index]*(f[IX(x,y)] - f[IX(x-1,y)]);          
      }
      else if (x ==0){
        adx = vx[index]*(f[IX(0,y)] - f[IX(N-1,y)]);  
      }
      else {  }
    
      if(y>0){
        ady = vy[index]*(f[IX(x,y)] - f[IX(x,y-1)]);          
      }
      else if (y ==0){
        ady = vy[index]*(f[IX(x,0)] - f[IX(x,M-1)]);  
      }
      else {  }
    
      ad.push_back(adx + ady);
    } //x
  }  //y

  return ad;
}


void Fluid::Den_evolve(){

/*
  float mass = 0;  
  for(int i=0;i<grid_size;i++){
      mass = mass + dens[i];
  }
  cout<<"Mass = "<<mass<<endl; 
*/
  
  vector<float> rhovx, rhovy, div_rhov; 

  for(int i=0;i<grid_size;i++){
    rhovx.push_back(dens[i]* vx[i] );
    rhovy.push_back(dens[i]* vy[i] );    
  }
//  div_rhov = Divergence(rhovx, rhovy, 1);
  
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



void Fluid::Vel_evolve(){
  
  vector<float>  vxvx, vxvy, vyvy, advection_vx, advection_vy;
  vector<float> ddx_P, ddy_P;
  vector<float>  laplace_vx, laplace_vy; 

  noslip();


  for(int i=0;i<grid_size;i++){
    vxvx.push_back(vx[i]*vx[i] );
    vxvy.push_back(vx[i]*vy[i] );
    vyvy.push_back(vy[i]*vy[i] );    
  }
  
  advection_vx = Divergence(vxvx, vxvy, 0);
  advection_vy = Divergence(vxvy, vyvy, 0);

  ddx_P = d_dx(P,0);   
  ddy_P = d_dy(P,0);      

  laplace_vx = laplacian(vx);
  laplace_vy = laplacian(vy);   

  float rhs;      
  for(int i=0; i<grid_size; i++){
  //- ddx_P[i]
  //- ddy_P[i]
  
    rhs  = -advection_vx[i] + (eta*laplace_vx[i] )/rho - ddx_P[i]/rho;
    temp_vx[i] = vx[i] + dt*rhs;

    rhs  = -advection_vy[i] + (eta*laplace_vy[i])/rho - ddy_P[i]/rho; ;
    temp_vy[i] = vy[i] + dt*rhs;  
    
   
    
  }
  vx = temp_vx;
  vy = temp_vy;


}



int it =0;

void Fluid::Evolve() {  
  
  vector<float> pressure_correction, dpc_dx, dpc_dy, div_v;
  float   error, p_avg;
  
  for(int iter=0; iter< 10; iter++){
    
    noslip();   
    Vel_evolve();
    noslip();
    pressure_correction = Divergence(vx, vy, 0);
    
    for(int i=0; i<grid_size; i++){
      pressure_correction[i] = -c_factor*pressure_correction[i];
    }
    dpc_dx = d_dx(pressure_correction, 0);
    dpc_dy = d_dy(pressure_correction, 0);
    
    for(int i=0; i<grid_size; i++){
      temp_vx[i] = vx[i] - dt_rho*dpc_dx[i]; 
      temp_vy[i] = vy[i] - dt_rho*dpc_dy[i]; 
      P[i] = P[i] + pressure_correction[i];
    }
    vx = temp_vx;
    vy = temp_vy;
    
    
    div_v = Divergence(vx, vy, 0 );
    
    error =0;

//    p_avg=0;
    for(int i=0; i<grid_size; i++){
      error = error + div_v[i]*div_v[i];
//      p_avg = p_avg + P[i];     
    }
    error = sqrt(error)/grid_size;
//    cout<<it<<"  avg_div_v : "<<error<<"  "<<p_avg<<endl;
    
    pressure_correction.clear();
    dpc_dx.clear(); 
    dpc_dy.clear();
    div_v.clear();
    
    if(error < 0.0001) break;
  } //pressure correction iteration

//if(it>1000) { Den_evolve();}
//   Den_evolve();
it++;
}


void Fluid::write_to_file(float t1, float t2){
  
  ofstream ff, bnd, ifd;
  
  ff.open("fluid_data_pipe10.txt");
  bnd.open("boundary_pipe.txt");
  ifd.open("instantanous_fluid_data_pipe10.txt");


  float t=t1, it=t2;
  int frame=0, index;
  
  for(int i=0; i<by.size(); i++ ){
    bnd<<setw(15)<<bx[i]<<setw(15)<<by[i]<<endl;    
  }  
  bnd.close();
  
  
  while (t <= t2){
  
    cout<<"t = "<<t<<endl;
    
      
    if(t>=t1 && t<=t2){
    
      for(int y=0; y<M; y++){
        for(int x=0; x<N; x++){
           index = IX(x,y);
           ff<<"    "<<frame<<"   "<<x<<"   "<<y<<"   "<<vx[index]<<"   "<<vy[index]<<"   "<<P[index]<<"    "<<t<<endl; 
        } //x
      }  //y

      if(t > it - dt && t< it + dt ){
            for(int y=0; y<M; y++){
        for(int x=0; x<N; x++){
           index = IX(x,y);
           ifd<<"    "<<frame<<"   "<<x<<"   "<<y<<"   "<<vx[index]<<"   "<<vy[index]<<"   "<<P[index]<<"    "<<t<<endl; 
        } //x
      }  //y


      
     }  //it
      

    } //t    
    
    Evolve();
    t = t + dt;
    frame++;
  }  //while  
  cout<<endl;

  ff.close();
  ifd.close();

}


















