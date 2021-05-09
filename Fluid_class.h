#ifndef FLUID_CLASS_H
#define FLUID_CLASS_H

#include <vector>

using namespace std;



class Fluid{

   public: 

    int N,M,grid_size;
    float eta = 0.5 , rho=1.0, Re;

    float del_x = 1;
    const float  dt= 0.01;

    float c_factor = rho*del_x*del_x/(4.0*dt), dt_rho = dt/rho;


    int IX(int x, int y);   


//dens is a scalar vector field over the screen 
    vector<float> temp_dens; 
    vector<float> dens; 
    
//    V the velocity field of the fluid 
    vector<float> temp_vx;
    vector<float> temp_vy;   
    vector<float> vx;
    vector<float> vy;
   
    vector<float> P;
    
//    w is the vorticitycity field of the fluid 
//    vector<float> temp_w;
 //   vector<float> w;
    
        
//  b is the set of boundary points
    vector<int> bx, by;
    
    Fluid();
    ~Fluid();
   
    void init_boundary();
    void init_density();
    void init_velocity();
 //   void init_vorticity();
    void noslip();
    void initialize();
 

    //void solve_steady_incompressible();
    vector<float> Advect(vector<float> f);
    void Den_evolve();
    void Vel_evolve();
    void Evolve();
    
    void write_to_file(float t1, float t2);
  
};  //  fluid class




#endif



/*    
// v = 60*v frame/sec
void Den_evolve(){
    
  int index1, index2;
  
  float x2,y2, temp;
      
  vector<float> tpx,tpy;    
      
  temp_dens = dens;  
 
   for(int i=0; i<px.size(); i++  ){
    index1 = IX( round(px[i]), round(py[i]) );  
    temp_dens[index1]=0; //  distruction 

    x2 = px[i] + Vx[index1]*dt;
    y2 = py[i] + Vy[index1]*dt;

    if(x2>=N-1 || x2<2){
      temp = -Vx[index1];
      x2 = px[i] + temp*dt;
    }

    if(y2>=N-1 || y2<2){
      temp = -Vy[index1];
      y2 = py[i] + temp*dt;
    }
    
    tpx.push_back(x2);
    tpy.push_back(y2);
    
    index2 = IX( round(x2), round(y2) ); 
    temp_dens[index2] = dens[index1];
    
  }

  dens =  temp_dens;
  px.clear();   py.clear();
  
  px = tpx;   py = tpy;
        
} 

*/









































