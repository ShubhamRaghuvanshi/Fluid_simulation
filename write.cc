#include"../Fluid_class.cc" 
#include "../Geometry.cc"
#include <unistd.h>

using namespace std;


vector<int> Bx, By;
float Vx[MG][NG]= {}, Vy[MG][NG]= {}, w[MG][NG]= {}, p[MG][NG]= {}; 
vector<float> VBX, VBY;


Fluid fluid;


int x1 = N/2-N/4, wy1 = M/2-M/4;
int x2 = N/2+N/4, wy2 = M/2-M/4;	
int x3 = N/2+N/4, wy3 = M/2+M/4;
int x4 = N/2-N/4, wy4 = M/2+M/4;


int main(int argc, char** argv)
{


	Line(x4, wy4, x3, wy3,	Bx, By);
//	Line(x4, wy4, x3-1, wy3,	topx, topy);

	int bs = Bx.size();
	for(int i=0; i<bs; i++){
		VBX.push_back(2);
		VBY.push_back(0);	
	} 

	Line(x1, wy1, x2, wy2,	Bx, By);
//	Line(x1, wy1, x2, wy2,	btmx, btmy);

	Line(x2, wy2, x3, wy3,	Bx, By);
//	Line(x2, wy2, x3, wy3,	rgtx, rgty);

	Line(x1, wy1, x4, wy4,	Bx, By);
//	Line(x1, wy1, x4, wy4,	lftx, lfty);

//	Line(x4-1, wy4+1, x3+1, wy3+1,	Bx, By);	
//	Line(x2+1, wy2-1, x3+1, wy3+1,	Bx, By);

//	Line(x1-1, wy1-1, x2+1, wy2-1,	Bx, By);
//	Line(x1-1, wy1, x4-1, wy4,	Bx, By);

	for(int i=0; i<Bx.size()-bs; i++){
		VBX.push_back(0);
		VBY.push_back(0);	
	} 

	for(int j=wy4+1; j<MG; j++){
		for(int i=0; i<NG; i++){
			Bx.push_back(i);
			By.push_back(j);
			VBX.push_back(0);
			VBY.push_back(0);	
		}
	}

	for(int j=0; j<wy1-1; j++){
		for(int i=0; i<NG; i++){
			Bx.push_back(i);
			By.push_back(j);
			VBX.push_back(0);
			VBY.push_back(0);	
		}
	}

	for(int j=wy2; j<=wy3; j++){
		for(int i=x2+1; i<NG; i++){
			Bx.push_back(i);
			By.push_back(j);
			VBX.push_back(0);
			VBY.push_back(0);	
		}
	}

	for(int j=wy1; j<=wy4; j++){
		for(int i=0; i<=x1-1; i++){
			Bx.push_back(i);
			By.push_back(j);
			VBX.push_back(0);
			VBY.push_back(0);	
		}
	}

//	for(int j=y4-3; j<=y4; j++){
//		for(int i=x1+1; i<x2-1;i ++){
//			fluid.svx[j][i]=1;
//			fluid.svy[j][i]=0;
//		}	
//	}

/*
	for(int j=wy4+1; j<MG; j++){
		for(int i=0; i<NG; i++){
			fluid.psibx.push_back(i);
			fluid.psiby.push_back(j);
			fluid.psib.push_back(0);
		}
	}

	for(int j=0; j<wy1-1; j++){
		for(int i=0; i<NG; i++){
			fluid.psibx.push_back(i);
			fluid.psiby.push_back(j);
			fluid.psib.push_back(0);
		}
	}

	for(int j=wy2; j<=wy3; j++){
		for(int i=x2+1; i<NG; i++){
			fluid.psibx.push_back(i);
			fluid.psiby.push_back(j);
			fluid.psib.push_back(0);
		}
	}

	for(int j=wy1; j<=wy4; j++){
		for(int i=0; i<=x1-1; i++){
			fluid.psibx.push_back(i);
			fluid.psiby.push_back(j);
			fluid.psib.push_back(0);		
		}
	}
*/

	fluid.initialize(Vx,Vy,p,Bx,By,VBX,VBY);

	char filename[25];
	sprintf(filename,"flow_%0.3f_%0.3f",fluid.Eu, fluid.Re);
	fluid.write_to_file(x1, wy1, x3-5, wy3-5, 0, 2.34, filename);



    return 0;
}































