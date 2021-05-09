#include <GL/glut.h>
#include"Fluid_class.h" 
#include"Geometry.h" 
#include <iostream>
#include <math.h>


using namespace std;



int main(){

  Fluid fluid;
  fluid.initialize();
  
  fluid.write_to_file(0,50);

  return 0;
}










