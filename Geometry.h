#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
using namespace std;

//Draws scalar field on 2D plane
void Draw_scalar_field(vector<float> f );
//mass of the fluid if f is density
float Mass(vector<float> f);
//Draws 2D vector field on 2D plane
void Draw_2dVector_Field( vector<float> fx, vector<float> fy);


//d/dx field with periodic boundary 
vector<float> d_dx(vector<float> f,  int schm );
//d/dy field with periodic boundary 
vector<float> d_dy(vector<float> f,  int schm );
//laplacian field with periodic boundary 
vector<float> laplacian(vector<float> f);


// divergence of vectr field f at point x y 
vector<float> Divergence(vector<float> fx, vector<float> fy, int schm);


//Draws white points
void Draw_points(vector<int> px, vector<int> py);
//Draws white points
void Draw_dens_points(vector<int> px, vector<int> py, vector<float> dens );

//Draws streamlines for the velocity field
void Draw_streamlines( vector<float> fx, vector<float> fy, float px, float py);


//Adds points
void point(int x , int y, vector<int> &bx, vector<int> &by );
//adds coordinates of pixels along a strainght line between two points
void Line( float x1, float y1, float x2, float y2, vector<int> &bx, vector<int> &by );
//points on a circle
void Circle( float c_x, float c_y, float r, vector<int> &bx, vector<int> &by, int n );
//points on rectangle
void rect( int x1, int y1, int x2, int y2, vector<int> &bx, vector<int> &by );
//points on solid rectangle
void solid_rect( int x1, int y1, int x2, int y2, vector<int> &bx, vector<int> &by );



#endif



















