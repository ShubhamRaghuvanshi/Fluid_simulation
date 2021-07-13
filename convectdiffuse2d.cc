
#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>

//#define float float 

#define px(i)  ( (i)  >= 0 ? (i)%N : N+(i) )
#define py(j)  ( (j)  >= 0 ? (j)%M : M+(j) )
 
using namespace std;

const int N=256, M=256;

float X1=-2, X2=2;
float Y1=-2, Y2=2;
float DX=(X2-X1)/(float(N)-1), DY = (Y2-Y1)/(float(M)-1); 

float RE=30000.0;
float D=1/RE;


float gauss_in(float A, float a, float x, float y, float x0, float y0){
	x = x-x0;
	y = y-y0;
	return A*exp(-a*(x*x+y*y)  );
}

float gauss_out(float A, float a, float x, float y, float t, float u, float v, float x0, float y0){
	x = x-x0;
	y = y-y0;
  float Ss = 1.0+4.0*a*D*t;	
	return 	( A/(1.0+4.0*a*D*t)*exp(-a*( (x-u*t)*(x-u*t) + (y-v*t)*(y-v*t) )/Ss ) );
}


//d/Dx field with periodic boundary using forward difference
void ddx_fd(float f[M][N], float DfDx[M][N]){
  for(int j=0; j<M; j++){
    for(int i=0; i<N-1; i++){
      	DfDx[j][i]=f[j][i+1]-f[j][i];
		}
	}		
  for(int j=0; j<M; j++){
    	DfDx[j][N-1]=f[j][0]-f[j][N-1];
	}	
} //fd


//d/Dy field with periodic boundary using forward difference
void ddy_fd(float f[M][N], float DfDy[M][N]){
  
  for(int j=0; j<M-1; j++){
    for(int i=0; i<N; i++){
      	DfDy[j][i]=f[j+1][i]-f[j][i];
		}
	}		
  for(int i=0; i<N; i++){
    	DfDy[M-1][i]=f[0][i]-f[M-1][i];
	}	
} //bd


//d/Dx field with periodic boundary using backward difference
void ddx_bd(float f[M][N], float DfDx[M][N]){
  for(int j=0; j<M; j++){
    for(int i=1; i<N; i++){
      	DfDx[j][i]=f[j][i]-f[j][i-1];
		}
	}		
  for(int j=0; j<M; j++){
    	DfDx[j][0]=f[j][0]-f[j][N-1];
	}	
} //bd


//d/Dy field with periodic boundary using backward difference
void ddy_bd(float f[M][N], float DfDy[M][N]){
  
  for(int j=1; j<M; j++){
    for(int i=0; i<N; i++){
      	DfDy[j][i]=f[j][i]-f[j-1][i];
		}
	}		
  for(int i=0; i<N; i++){
    	DfDy[0][i]=f[0][i]-f[M-1][i];
	}	
} //bd


//d/Dx field with periodic boundary using central difference
void ddx_cd(float f[M][N], float DfDx[M][N]){

  for(int j=0; j<M; j++){
    for(int i=1; i<N-1; i++){
      	DfDx[j][i]=f[j][i+1]-f[j][i-1] ;
		}
	}	
  for(int j=0; j<M; j++){
 		DfDx[j][0] = f[j][1] - f[j][N-1];  
		DfDx[j][N-1] = f[j][0] - f[j][N-2];  
	}	
} //cd


//d/Dy field with periodic boundary using central difference
void ddy_cd(float f[M][N], float DfDy[M][N]){
  
  for(int j=1; j<M-1; j++){
    for(int i=0; i<N; i++){
      	DfDy[j][i]=f[j+1][i]-f[j-1][i];
		}
	}	
  for(int i=0; i<N; i++){
  	DfDy[0][i]=f[1][i]-f[M-1][i] ;
  	DfDy[M-1][i]=f[0][i]-f[M-2][i] ;
	}	
} //cd

void ddf_dxx(float f[M][N], float g[M][N]){
  for(int j=0; j<M; j++){
    for(int i=1; i<N-1; i++){
			g[j][i] =  f[j][i+1] - 2.0*f[j][i] + f[j][i-1] ; 
    }
	}
  for(int j=0; j<M; j++){
			g[j][0]   = f[j][1] - 2.0*f[j][0] + f[j][N-1] ; 
			g[j][N-1] = f[j][0] - 2.0*f[j][N-1] + f[j][N-2] ; 
	}
}

void ddf_dyy(float f[M][N], float g[M][N]){
  for(int j=1; j<M-1; j++){
    for(int i=0; i<N; i++){
      g[j][i] =  f[j+1][i] - 2.0*f[j][i] + f[j-1][i] ;
    }
	}
  for(int i=0; i<N; i++){
    g[0][i] =  f[1][i] - 2.0*f[0][i] + f[M-1][i] ;
    g[M-1][i] = f[0][i] - 2.0*f[M-1][i] + f[M-2][i] ;
  }
}

void Advect(float f[M][N], float vx[M][N], float vy[M][N], float S[M][N], float Dt ){

	float dfdx[M][N], dfdy[M][N]; 
	ddx_bd(f,dfdx);
	ddy_bd(f,dfdy);

	float alpx = Dt/DX, alpy = Dt/DY;// beta = D*Dt/(DX*DX) 
  for(int j=0; j<M; j++){
    for(int i=0; i<N; i++){
			f[j][i] = f[j][i] - alpx*vx[j][i]*dfdx[j][i] - alpy*vy[j][i]*dfdy[j][i] - Dt*S[j][i];
		}
	}		
}

float Max(float f[M][N]){
	float U=0;
	for(int j=0; j<M; j++){
	  for(int i=0; i<N; i++){
			if(U < abs(f[j][i]) ) { U = abs(f[j][i]);  }
		}
	}
	return U;
}



float Max(float fx[M][N], float fy[M][N]){
	float U=0,u=0;
	for(int j=0; j<M; j++){
	  for(int i=0; i<N; i++){
			u = fx[j][i]*fx[j][i]	+ fy[j][i]*fy[j][i];
			if(U<u) {U=u;}
		}
	}
	return sqrt(U);
}

bool chk=false;
void check(float u[M][N], float v[M][N], float Dt){
	if(!chk){
		float U=Max(u,v);
		if( Dt < 2.0*D/(U*U) && DX > 2.0*D*Dt/U ){
			chk = true;
		}
		else{
		  cout<<"Warning : Solution might be unstable"<<endl;
		  cout<<"dx = "<<DX<<" ,   dt = "<<Dt<<" ,    U = "<<U<<endl;
			cout<<"For stable solution use:  dt < "<<2.0*D/(U*U);
			cout<<",  dx >"<<2.0*D*Dt/U<<endl<<endl;  	  
			chk = true;
	  }
	}	  
}


void DeltaF(float k[M][N], float f[M][N], float Vx[M][N], float Vy[M][N], float S[M][N], float delt){

	float Ax =delt/(2.0*DX), Ay =delt/(2.0*DY), Bx = D*delt/(DX*DX), By = D*delt/(DY*DY);

	float dfdx[M][N], dfdy[M][N], ddfdxx[M][N], ddfdyy[M][N];

	ddx_cd(f,dfdx);
	ddy_cd(f,dfdy);
	ddf_dxx(f, ddfdxx);
	ddf_dyy(f, ddfdyy);

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			k[j][i] =  -Ax*Vx[j][i]*dfdx[j][i] -Ay*Vy[j][i]*dfdy[j][i] +Bx*ddfdxx[j][i] +By*ddfdyy[j][i] -delt*S[j][i] ;
		}
	}			
}


void add( float a[M][N], float b[M][N] ){
	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			a[j][i] = a[j][i] + b[j][i];  
		}
	}			
}

void set_bnd(float f0[M][N], vector<int> bx, vector<int> by, vector<float> F){
	if(bx.size() == by.size() && bx.size() == F.size() ){
		for(int i=0; i<bx.size(); i++){
			f0[ by[i] ][ bx[i] ] = F[i];
		}
	}
	else {
		cout<<"Array size dont match : "<<bx.size()<<"	"<<by.size()<<"	"<<F.size()<<endl;	
	}
}


void Explicit_FTCS(float f[M][N], float Vx[M][N], float Vy[M][N], float S[M][N], vector<int> bx, vector<int> by, vector<float> fb,  float Dt){

	float k[M][N];

	check(Vx, Vy, Dt);
	set_bnd(f, bx, by, fb);
	DeltaF(k, f, Vx, Vy, S, Dt);
	add(f, k);
	set_bnd(f, bx, by, fb);
}
	

void Explicit_FTCS_rk2(float f[M][N], float Vx[M][N], float Vy[M][N], float S[M][N], vector<int> bx, vector<int> by, vector<float> fb, float Dt){

	float k1[M][N], k2[M][N], f1[M][N];

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f1[j][i] = f[j][i];
		}
	}

	check(Vx, Vy, Dt);
	set_bnd(f1, bx, by, fb);
	DeltaF(k1, f1, Vx, Vy, S, Dt);
	add(f1, k1);
	set_bnd(f1, bx, by, fb);
	DeltaF(k2, f1, Vx, Vy, S, Dt);

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f[j][i] = f[j][i] + (k1[j][i]+k2[j][i])/2.0;
		}
	}
	set_bnd(f, bx, by, fb);
}

void BTCS(float f[M][N], float Vx[M][N], float Vy[M][N], float guess[M][N], float S[M][N], vector<int> bx, vector<int> by, vector<float> fb, float Dt ){

	float alpx =Dt/(2.0*DX), alpy =Dt/(2.0*DY);
	float Ax, Ay,  Bx = D*Dt/(DX*DX), By = D*Dt/(DY*DY);
	float err, f0[M][N], f1[M][N], term;

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f0[j][i] = guess[j][i];
		}
	}
	set_bnd(f0, bx, by, fb);


	for( int itr=0; itr<1000; itr++ ){
		err =0;		
		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				Ax = Vx[j][i]*alpx; 	Ay = Vy[j][i]*alpy;
				term =  f[j][i] -(Ax-Bx)*f0[j][px(i+1)] -(Ay-By)*f0[py(j+1)][i] +(Ax+Bx)*f0[j][px(i-1)] +(Ay+By)*f0[py(j-1)][i] -Dt*S[j][i]; 
				f1[j][i] = term/(1.0+2.0*Bx+2.0*By);
			}
		}
		set_bnd(f1, bx, by, fb);

		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				err = err + abs( f1[j][i]-f0[j][i]);				
			}
		}
		err = err/(N*M);
		if(err < 0.000001) { 
		//cout<<"btcs itr : "<<itr<<endl;	
		break;
		}

		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				f0[j][i]=f1[j][i];
			}
		}
	}	

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f[j][i]=f1[j][i];
		}
	}
}

void predictor_corrector(float f[M][N], float Vx[M][N], float Vy[M][N], float S[M][N], vector<int> bx, vector<int> by, vector<float> fb, float dt ){

	float guess[M][N], k[M][N] ;

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			guess[j][i]=f[j][i];
		}
	}
//  set_bnd(guess, bx, by, fb);
//	DeltaF(k, guess, Vx, Vy, S, dt);
//	add(guess, k);
//  set_bnd(guess, bx, by, fb);


  Explicit_FTCS_rk2(guess , Vx, Vy, S, bx, by, fb, dt);
	BTCS(f, Vx, Vy, guess, S,  bx, by, fb, dt );
}

void explicit_implicit(float f[M][N], float Vx[M][N], float Vy[M][N], float S[M][N], vector<int> bx, vector<int> by, vector<float> fb, float delt ){

	float k[M][N], fh[M][N], guess[M][N] ;

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			fh[j][i]=f[j][i];
			guess[j][i]=f[j][i];
		}
	}
	set_bnd(fh, bx, by, fb);
	DeltaF(k, fh, Vx, Vy, S, delt/2.0);
	add(fh, k);
	set_bnd(fh, bx, by, fb);

//	DeltaF(k, guess, Vx, Vy, S, delt);
//	add(guess, k);
	Explicit_FTCS_rk2(guess, Vx, Vy, S, bx, by, fb, delt);

	BTCS(fh, Vx, Vy, guess, S, bx, by, fb, delt/2.0 );

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f[j][i] = fh[j][i];
		}
	}
}


//D^2 f = -g 
void poisson_solve(float f[M][N], float g[M][N], vector<int> bx, vector<int> by, vector<float> fb ){

	float f0[M][N], f1[M][N], DXX= DX*DX, DYY = DY*DY;
	float A = DYY*DXX/( 2.0*(DYY + DXX) ), B = DYY/( 2.0*(DYY + DXX) ),C = DXX/( 2.0*(DYY + DXX) ); 
	float error;


	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f0[j][i] = f[j][i];
		}
	}
	
	set_bnd(f0, bx, by , fb);		

	for(int itr=0; itr<1000; itr++){
		error =0;
		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				f1[j][i] = g[j][i]*A + (f0[j][px(i+1)] + f0[j][px(i-1)])*B+ (f0[py(j+1)][i] + f0[py(j-1)][i])*C;
			}	
		}
		set_bnd(f1, bx, by , fb);	
		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				error = error + abs( f1[j][i] - f0[j][i] ); 
			}	
		}
		error = error/(M*N);
	//	cout<<error<<endl;
		if(error < 0.00001) { cout<<"poisson itr : "<<itr<<endl; break; }
		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				f0[j][i] = f1[j][i];
			}
		}
	}
//	set_bnd(f1, bx, vy , F);	
	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			f[j][i] = f1[j][i];
		}
	}

}











