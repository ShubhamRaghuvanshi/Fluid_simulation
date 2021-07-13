#include<iostream>
#include <fstream>
#include"convectdiffuse1d.cc"

using namespace std;

int main(){

	dtype x=x1, t, phi[N], dphi[N], x0=0.5;
	dtype A=1, a=200, t1=0, t2=0.3, u=1	;
	int  ti=0, tN = (t2-t1)/DT;
	dtype U[N];
	for(int i=0; i<N; i++){
		U[i]=u;
	}
		
	char filename[30];
	
	sprintf(filename, "phi_%d.txt", ti );
  ofstream output;
  output.open(filename);
	for(int i=0; i<N; i++){
		phi[i] = gauss_in(A, a, x, x0);
		output<<ti<<"		"<<x<<"		"<<phi[i]<<endl;
		x = x1 + (i+1)*DX;
	}
	output.close();
	
	for(ti=0; ti<tN; ti++){
//	Explicit_FTCS(phi , U, DT, dphi);
//	Explicit_FTCS_lax(phi , U, DT, dphi);
// 	Explicit_FTCS_rk2(phi, U, DT);	
 //	Explicit_FTCS_rk4(phi, U, DT);	
	//		Upwind_convection(phi , U, DT);
//	Implicit_FTCS(phi , U, phi, DT);
//	Predictor_corrector(phi, U, DT);
//	Explicit_Implicit(phi, U, DT);

		t = t1 + (ti+1)*DT;	
	}

	sprintf(filename, "phi_%d.txt", ti );
  output.open(filename);
  x=x1;
  dtype error=0, analsol, I1=0, I2=0;
	for(int i=0; i<N; i++){
//		output<<ti<<"		"<<x<<"		"<<gauss_out(A, a, x, 0.3, u, x0)<<"		"<<gauss_out(A, a, x, 0.3, u, x0)<<endl;
		analsol = gauss_out(A, a, x, t2, u, x0);
		output<<ti<<"		"<<x<<"		"<<analsol<<"		"<<phi[i]<<endl;
		if(abs(analsol)>1e-4)
			error = error + abs( (phi[i] -  analsol)/analsol );
		I1 = I1 + analsol;
		I2 = I2 + phi[i];
		x = x1 + (i+1)*DX;
	}
	output.close();
	error = error/N;
	cout<<"int : "<<I1<<"		"<<I2<<endl;

  FILE *gnuplot;
  gnuplot = popen("gnuplot -persist", "w");
  if(gnuplot != NULL){
      //Uncomment next two lines to get image of the output
  //   fprintf(gnuplot, "set terminal png size 640,480\n" );
   //   fprintf(gnuplot, "set output 'i7.png' \n");
			fprintf(gnuplot, "set title 'A=%0.2f, a=%0.2f, D=%0.2f, u=%0.2f, dt=%0.3f' \n", A, a, D, u, DT );
      fprintf(gnuplot, "set xlabel 'x with dx = %0.3f' \n", DX );      
    fprintf(gnuplot, "plot 'phi_%d.txt' using 2:3 w l title 'phi(x,%0.3f)', 'phi_%d.txt' using 2:4 w l title 'Explicit-Implicit'\n", ti, t, ti);        
 
  }
//  cout<<"Created image file : p5_"<<itr<<".png "<<endl;



	return 0;
}
