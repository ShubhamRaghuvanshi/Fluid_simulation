#ifndef FLUID_CLASS_H
#define FLUID_CLASS_H


#include "convectdiffuse2d.cc"

class Fluid{

	public: 

		float nu = 10 , rho=1.0, Re=RE, Eu=1;

		float t=0, dt = 0.005;
	
		float v_factor = DX*DY/(4.0*(DX+DY));
		float p_factor = v_factor*rho/dt;

		float vx[M][N]   = {};
		float vy[M][N]   = {};      
		float P[M][N]    = {};

		float w[M][N]={}, psi[M][N]={}; 
		float sw[M][N]={}, svx[M][N], svy[M][N]; 


		//  (bx,by) is the set of boundary points
		vector<int> bx, by, psibx, psiby, nbx, nby;
		vector<float> vxb, vyb, psib, wb;

/*
		// velocity source (if any)
		vector<int> xvs, yvs;
		vector<float> vsx, vsy;

		// density source (if any)
		vector<int> xds, yds;
		vector<float> ds;
		float dens[M][N] = {};
*/


		Fluid();
		~Fluid();

		void init_boundary(vector<int> Bx, vector<int> By, vector<float> VXb	, vector<float> VYb);
		
		void init_p(float p[M][N]);
		void init_v(float fx[M][N], float fy[M][N]);
		void init_w();
		void init_psi(float Psi[M][N]);

		void vort();

		void initialize(float fx[M][N], float fy[M][N], float p[M][N], vector<int> Bx, vector<int> By, vector<float> vBx, vector<float> vBy );

		void setbnd_v();
		void setbnd_w();
		void setbnd_psi();

	

		float func_wbnd(int n, int m, int BX, int BY );
		void vorticity_stream();
		void pressure_corrector();


		void write_to_file(int x1, int y1, int x2, int y2, float t1, float t2, char filename[]);
};  //  fluid class


void read_from_file(float &T, float fx[M][N], float fy[M][N], float p[M][N], vector<int> &Bx, vector<int> &By, char filename[]);
#endif







































