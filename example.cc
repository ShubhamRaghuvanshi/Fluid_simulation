#include<iostream>
#include <fstream>
#include"convectdiffuse2d.cc"

using namespace std;




int main(){

	datatype x=X1, y=Y1, t, phi[M][N], x0=0.5, y0=0.5;
	datatype A=1, a=200, t1=0, t2=0.3, u=1, v=0	;
	
	int  ti=0, tN = (t2-t1)/DT;
	datatype vx[M][N], vy[M][N];

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			vx[j][i]=u;
			vy[j][i]=v;
		}
	}
			
	char filename[30];
	sprintf(filename, "phi_%d.txt", ti );


  ofstream output;

/*
  output.open(filename);

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			phi[j][i] = gauss_in(A, a, x, y, x0, y0);
			output<<ti<<"		"<<x<<"	"<<y<<"		"<<phi[j][i]<<endl;
			x = X1 + (i+1)*DX;
		}
		y = Y1 + (j+1)*DX;
	}
	output.close();
*/	
//	for(ti=0; ti<tN; ti++){
//	Explicit_FTCS(phi, vx, vy, DT);
//		Explicit_FTCS_rk2(phi, vx, vy, DT);
	//	BTCS(phi, vx, vy, phi, DT);
//	predictor_corrector(phi, vx, vy, DT);	
//	explicit_implicit(phi, vx, vy, DT);	
//		t = t1 + (ti+1)*DT;	
//	}


	sprintf(filename, "phi_%d.txt", ti );
  output.open(filename);
  x=X1;
	y =Y1;
  datatype error=0, analsol, I1=0, I2=0;
	for(int j=0; j<M; j++){

		for(int i=0; i<N; i++){
//		output<<ti<<"		"<<x<<"		"<<gauss_out(A, a, x, 0.3, u, x0)<<"		"<<gauss_out(A, a, x, 0.3, u, x0)<<endl;
			
			output<<ti<<"		"<<x<<"		"<<y<<"		"<<gauss_out(A, a, x,y, 0, u, v, x0, y0)<<"		"<<gauss_out(A, a, x,y, t2, u, v, x0, y0)<<endl;
//		if(abs(analsol)>1e-4)
			error = error + abs( phi[j][i] -  analsol);
//		I1 = I1 + analsol;
//		I2 = I2 + phi[i];
			x = X1 + (i+1)*DX;
		}
		y = Y1 + (j+1)*DY;
	}
	

	output.close();
	error = error/(N*M);
//	cout<<"int : "<<I1<<"		"<<I2<<endl;


  FILE *gnuplot;
  gnuplot = popen("gnuplot -persist", "w");
  if(gnuplot != NULL){
      //Uncomment next two lines to get image of the output
  //   fprintf(gnuplot, "set terminal png size 640,480\n" );
    //  fprintf(gnuplot, "set output 'I7.png' \n");
			fprintf(gnuplot, "set title 'A=%0.2f, a=%0.2f, D=%0.2f, u= %0.2f, v=%0.2f ' \n", A, a, D, u, v );
//			fprintf(gnuplot, "set title 'Explicit implicit, error = %0.5f,  t = %0.4f '\n", error, t2 );
      fprintf(gnuplot, "set xlabel 'x with dx = %0.3f' \n", DX ); 
      fprintf(gnuplot, "set ylabel 'y with dy = %0.3f' \n", DY );
			fprintf(gnuplot, "set view map\n");
	    fprintf(gnuplot, "splot 'phi_%d.txt' using 2:3:4 w p pointtype 5 palette title '', 'phi_%d.txt' using 2:3:5 w p pointtype 5 palette title ''\n",ti, ti);        
 	//    fprintf(gnuplot, "splot 'phi_%d.txt' using 2:3:4 w l title ''\n", ti);        
 
 
//    fprintf(gnuplot, "plot 'phi_%d.txt' using 2:3 w l title 'phi(x,%0.3f)', 'phi_%d.txt' using 2:4 w l title 'Explicit-Implicit'\n", ti, t, ti);        
 
  }
//  cout<<"Created image file : p5_"<<itr<<".png "<<endl;



	return 0;
}
