#include <iostream>
#include <math.h>
#include <GL/glut.h>

#define dtype float
 
using namespace std;

const int N=100, M=128, SCALE=2;

dtype x1=-1, x2=3;
dtype Y1=0, Y2=0.25;
dtype DX=(x2-x1)/(float(N)-1), D=0, DT=0.01;


dtype gauss_in(dtype A, dtype a, dtype x, dtype x0){
	x = x-x0;
	return A*exp(-a*x*x);
}

dtype gauss_out(dtype A, dtype a, dtype x, dtype t, dtype u, dtype x0){
	x = x-x0;
  dtype Ss = 1.0+4.0*a*D*t;	
	return 	(A/sqrt(1.0+4.0*a*D*t))*exp(-a*(x-u*t)*(x-u*t)/Ss );
}

dtype impulse_resp(dtype a, dtype x, dtype t, dtype u, dtype x0){
	x = x- x0;
	return 	(1.0/sqrt(4.0*M_PI*D*t))*exp(-(x-u*t)*(x-u*t)/(4.0*D*t) );
}


void Advect_I(dtype f[N], dtype vx[N]){
  for(int i=1; i<N; i++){
		f[i] = f[i] - vx[i]*DT/DX*(f[i] - f[i-1]) ;
	}

	f[0] = f[0] - vx[0]*DT/DX*(f[0] - f[N-1]) ;

}

void Advect_C(dtype f[N], dtype vx[N],  dtype DT ){
  for(int i=1; i<N; i++){
		f[i] = f[i] - DT/DX*( vx[i]*f[i] - vx[i-1]*f[i-1]) ;
	}
	f[0] = f[0] - DT/DX*(vx[0]*f[0] - vx[N-1]*f[N-1]) ;
}


bool chk=false;
void check(dtype u[N]){
	if(!chk){
		dtype U=0;
		for(int i=0; i<N; i++){
			if(U < abs(u[i])) {U = abs(u[i]); }
		}

		if( DT < 2.0*D/(U*U) && DX > U*DT ){
			chk = true;
		}
		else{
		  cout<<"Warning : Solution might be unstable"<<endl;
		  cout<<"dx = "<<DX<<" ,   dt = "<<DT<<" ,    U = "<<U<<endl;
			cout<<"For stable solution use:  dt < "<<2.0*D/(U*U);
			cout<<",  dx >"<<U*DT<<endl<<endl;  	  
			chk = true;
	  }
	}	  
}


void DeltaF(dtype k[N], dtype f[N], dtype Vx[N], dtype delt){

	dtype alp =delt/(2.0*DX),  beta = D*delt/(DX*DX);

	for(int i=1; i<N-1; i++){		
		k[i] = - alp*Vx[i]*( f[i+1]- f[i-1] )  + beta*( f[i+1] - 2.0*f[i] + f[i-1]);
	}
	k[0] = 	- alp*Vx[0]   *( f[1]- f[N-1] )  + beta*( f[1] - 2.0*f[0]    + f[N-1]);
	k[N-1] = - alp*Vx[N-1]*( f[0]- f[N-2] )  + beta*( f[0] - 2.0*f[N-1] + f[N-2]);
}


/*
void DeltaF(dtype k[N], dtype f[N], dtype Vx[N], dtype delt){

	dtype alp =delt/(12.0*DX),  beta = D*delt/(12.0*DX*DX);

	for(int i=2; i<N-2; i++){		
		k[i] = - alp*Vx[i]*( -f[i+2]+ 8.0*f[i+1]- 8.0*f[i-1] + f[i-2])  + beta*(  -f[i+2]+ 16.0*f[i+1]  -30.0*f[i] +16.0*f[i-1] - f[i-2] ) ;
	}
	k[0] = - alp*Vx[0]*( -f[2]+ 8.0*f[1]- 8.0*f[N-1] + f[N-2])  + beta*(  -f[2]+ 16.0*f[1]  -30.0*f[0] +16.0*f[N-1] - f[N-2] ) ;
	k[1] = - alp*Vx[1]*( -f[2]+ 8.0*f[2]- 8.0*f[0] + f[N-1])  + beta*(  -f[3]+ 16.0*f[2]  -30.0*f[1] +16.0*f[0] - f[N-1] ) ;

	k[N-2] = - alp*Vx[N-2]*( -f[0]+ 8.0*f[N-1]- 8.0*f[N-3] + f[N-4])  + beta*(  -f[0]+ 16.0*f[N-1]  -30.0*f[N-2] +16.0*f[N-3] - f[N-4] ) ;
	k[N-1] = - alp*Vx[N-1]*( -f[1]+ 8.0*f[0]- 8.0*f[N-2] + f[N-3])  + beta*(  -f[1]+ 16.0*f[0]  -30.0*f[N-1] +16.0*f[N-2] - f[N-3] ) ;
}
*/


void add(dtype a[N], dtype b[N]){
	for(int i=0; i<N; i++){		
		a[i] = a[i] + b[i]; 	
	}	
}

void Explicit_FTCS(dtype f[N], dtype Vx[N], dtype Dt, dtype k[N]){
	check(Vx);
	DeltaF(k, f, Vx, Dt);
	add(f, k);
}
	
void Upwind_convection(dtype f[N], dtype Vx[N], dtype Dt){

	dtype alp =Dt/DX,  beta = D*Dt/(DX*DX);
	dtype k[N];
	check(Vx);

	for(int i=1; i<N-1; i++){		
		k[i] = -alp*Vx[i]*( f[i]- f[i-1] )  + beta*( f[i+1] - 2.0*f[i] + f[i-1]);
	}
	k[0] = 	-alp*Vx[0]*( f[0]- f[N-1] )  + beta*( f[1] - 2.0*f[0]    + f[N-1]);
	k[N-1] = -alp*Vx[N-1]*( f[N-1]- f[N-2] )  + beta*( f[0] - 2.0*f[N-1] + f[N-2]);

	add(f, k);	
}



void Explicit_FTCS_lax(dtype f[N], dtype Vx[N], dtype Dt,	dtype k[N]){

	dtype avg_f[N];
	for(int i=1; i<N-1; i++){		
		avg_f[i] =  (f[i+1]+f[i-1])/2.0;
	}
	avg_f[0] = 	(f[1]+f[N-1])/2.0; 
	avg_f[N-1] = (f[0]+f[N-2])/2.0;


	check(Vx);
	DeltaF(k, f, Vx, Dt);
	for(int i=0; i<N; i++){		
		f[i] = avg_f[i];
	}
	add(f, k);
}

	
void Explicit_FTCS_rk2(dtype f[N], dtype Vx[N], dtype Dt){

	dtype k1[N], k2[N], f1[N];
	check(Vx);

	for(int i=0; i<N; i++){		
		f1[i] = f[i]; 	
	}	
	DeltaF(k1, f, Vx, Dt);
	add(f1, k1);
	DeltaF(k2, f1, Vx, Dt);

	for(int i=0; i<N; i++){
		f[i] = f[i] + ( k1[i]+k2[i])/2.0;
	}
}

void Explicit_FTCS_rk4(dtype f[N], dtype Vx[N], dtype Dt){

	dtype k1[N], k2[N], k3[N], k4[N]; 
	dtype f1[N], f2[N], f3[N];
	check(Vx);

	for(int i=0; i<N; i++){		
		f1[i] = f[i];
		f2[i] = f[i]; 	
		f3[i] = f[i]; 	
	}	

	DeltaF(k1, f, Vx, Dt/2.0); // k1 = f(yn)
	add(f1, k1);            // f1 = yn +hk1/2;

	DeltaF(k2, f1, Vx, Dt/2.0); // k2 = f(yn +hk1/2)
	add(f2, k2);            // f2 = yn +hk2/2;

	DeltaF(k3, f2, Vx, Dt); // k3 = f(yn +hk2/2)
	add(f3, k3);            // f3 = yn +hk3;

	DeltaF(k4, f3, Vx, Dt); //k4 = f(yn +hk3)

	for(int i=0; i<N; i++){
		f[i] = f[i] + ( 2.0*k1[i]+4.0*k2[i]+2.0*k3[i]+k4[i])/6.0;
	}
}


void Implicit_FTCS(dtype f[N], dtype Vx[N], dtype guess[N], dtype Dt){

	dtype alp =Dt/(2.0*DX),  beta = D*Dt/(DX*DX);
	dtype alpha, err, f0[N], f1[N];

	for(int i=0; i<N; i++){		
		f0[i] = guess[i];
	}

	for( int itr=0; itr<10000; itr++ ){

		for(int i=1; i<N-1; i++){		
			alpha = Vx[i]*alp;
			f1[i] = (f[i] - (alpha-beta)*f0[i+1] +  (alpha+beta)*f0[i-1])/(1.0+2.0*beta);
		}		
		alpha = Vx[0]*alp;
		f1[0] = (f[0] - (alpha-beta)*f0[1] +  (alpha+beta)*f0[N-1])/(1.0+2.0*beta);
		alpha = Vx[N-1]*alp;
		f1[N-1] = (f[N-1] - (alpha-beta)*f0[0] +  (alpha+beta)*f0[N-2])/(1.0+2.0*beta);

		err =0;		
		for(int i=0; i<N; i++){		
			err = err + abs ( f0[i] -f1[i]  );
		}		
		err = err/N;
		if(err < 0.00001) { break;}

		for(int i=0; i<N; i++){		
			f0[i] = f1[i];
		}
	}	

	for(int i=0; i<N; i++){	
		f[i] = f1[i];
	}		
	
}


void Predictor_corrector(dtype f[N], dtype Vx[N], dtype Dt){

	dtype guess[N], k[N];
	for(int i=0; i<N; i++){		
		guess[i] = f[i];
	}
	DeltaF(k, guess, Vx, Dt);
	add(guess, k);	
	Implicit_FTCS(f,Vx, guess, Dt);	
}


void Explicit_Implicit(dtype f[N], dtype Vx[N], dtype Dt){

	dtype k[N], fh[N], guess[N];

	for(int i=0; i<N; i++){		
		fh[i] = f[i];
		guess[i] = f[i];
	}

	DeltaF(k, fh, Vx, Dt/2.0);
	add(fh, k);	

	DeltaF(k, guess, Vx, Dt);
	add(guess, k);	

	Implicit_FTCS(fh, Vx, guess, Dt/2.0);

	for(int i=0; i<N; i++){		
		f[i]= fh[i];
	}
}















