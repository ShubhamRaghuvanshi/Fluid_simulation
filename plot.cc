#include<iostream>
#include <fstream>

using namespace std;


int main(){

	char filename[25];
	float EU=1, RE=400;
	sprintf(filename,"flow_%0.3f_%0.3f", EU, RE);

  FILE *gnuplot;
  gnuplot = popen("gnuplot -persist", "w");
  if(gnuplot != NULL){
      //Uncomment next two lines to get image of the output
     fprintf(gnuplot, "set terminal png size 640,480\n" );
      fprintf(gnuplot, "set output 'w.png' \n");
	//		fprintf(gnuplot, "set title 'A=%0.2f, a=%0.2f, D=%0.2f, u= %0.2f, v=%0.2f ' \n", A, a, D, u, v );
//			fprintf(gnuplot, "set title 'Explicit implicit, error = %0.5f,  t = %0.4f '\n", error, t2 );
      fprintf(gnuplot, "unset xtics \n" ); 
      fprintf(gnuplot, "unset ytics\n" );
			fprintf(gnuplot, "set view map\n");
			fprintf(gnuplot, "set title 'Vorticity'\n" );
	    fprintf(gnuplot, "splot 'flow_1.000_400.000_final.txt'  using 3:4:8 w p pointtype 5 palette title ''\n");    
			fprintf(gnuplot, "set title 'Stream function'\n" );
	    fprintf(gnuplot, "set output 'psi.png' \n");
    	fprintf(gnuplot, "splot 'flow_1.000_400.000_final.txt'  using 3:4:9 w p pointtype 5 palette title ''\n");    
    
 	//    fprintf(gnuplot, "splot 'phi_%d.txt' using 2:3:4 w l title ''\n", ti);        
 	
 
 
//    fprintf(gnuplot, "plot 'phi_%d.txt' using 2:3 w l title 'phi(x,%0.3f)', 'phi_%d.txt' using 2:4 w l title 'Explicit-Implicit'\n", ti, t, ti);        
 
  }
//  cout<<"Created image file : p5_"<<itr<<".png "<<endl;



	return 0;
}
