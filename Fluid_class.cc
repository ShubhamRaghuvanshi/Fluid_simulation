#include"Fluid_class.h"
#include<cstring>
#include<fstream>

using namespace std;




Fluid::Fluid(){
}

Fluid::~Fluid(){}; //distructor


void Fluid::init_boundary(vector<int> Bx, vector<int> By, vector<float> VXb	, vector<float> VYb){

	if(Bx.size() == By.size()){
		cout<<"Init boundary"<<endl;
		for(int i=0; i<Bx.size(); i++){
			bx.push_back(Bx[i]);
			by.push_back(By[i]);

			vxb.push_back(VXb[i]);
			vyb.push_back(VYb[i]);

			psibx.push_back(Bx[i]);
			psiby.push_back(By[i]);

			psib.push_back(0);
		}
	}
	
	else{
		cout<<"Boundary point error : Array not same size : "<<Bx.size()<<"	"<<By.size()<<endl;
	} 
//	not_bnd(nbx, nby);
}	

void Fluid::init_p(float p[M][N]){

	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			P[j][i] = p[j][i];
//			P[j][i]= gauss_out(1, 200, X1+(i+1)*DX, Y1+(j+1)*DY, 0, 1, 0, 0, 0);
		}		
	}
}


void Fluid::init_v(float fx[M][N], float fy[M][N]){
	
	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			vx[j][i] = fx[j][i];
			vy[j][i] = fy[j][i];
		}		
	}
}


void Fluid::init_w(){
	
	float dxvy[M][N], dyvx[M][N];	
	ddx_cd( vy, dxvy );
  ddy_cd( vx, dyvx );

  for(int j=0; j<M; j++){
    for(int i=0; i<N; i++){
    	w[j][i] = dxvy[j][i]/(2.0*DX) - dyvx[j][i]/(2.0*DY);
		}
	}
}

void Fluid::vort(){
	float dxvy[M][N], dyvx[M][N];	
	ddx_cd( vy, dxvy );
  ddy_cd( vx, dyvx );

  for(int j=0; j<M; j++){
    for(int i=0; i<N; i++){
    	w[j][i] = dxvy[j][i]/(2.0*DX) - dyvx[j][i]/(2.0*DY);
		}
	}
}


void Fluid::init_psi(float Psi[M][N]){

  for(int j=0; j<M; j++){
    for(int i=0; i<N; i++){
    	psi[j][i] = Psi[j][i]; 
		}
	}
}


float Fluid::func_wbnd(int n, int m, int BX, int BY ){

	 return ( 2.0*(psi[BY][BX]- psi[BY+m][BX+n])/(DX*DX) + 2.0*( m*vx[BY][BX] - n*vy[BY][BX] )/DX );   
}


vector <int> topx, topy, btmx, btmy, rgtx, rgty, lftx, lfty;
// applies noslip on boundary
void Fluid::setbnd_v(){
  for(int i=0; i<bx.size(); i++ ){
    vx[by[i]][bx[i]] = vxb[i];
    vy[by[i]][bx[i]] = vyb[i];
   }  
}

void Fluid::setbnd_w(){
	for(int i=0; i<topx.size(); i++){
		 	 w[topy[i]][topx[i]] = func_wbnd(0, -1, topx[i], topy[i] );	 	
	 }	

   for(int i=0; i<btmx.size(); i++){
		 	 w[btmy[i]][btmx[i]] = func_wbnd(0, 1, btmx[i], btmy[i] );	 	
	 }	

	 for(int i=0; i<rgtx.size(); i++){
		 	 w[rgty[i]][rgtx[i]] = func_wbnd(-1, 0, rgtx[i], rgty[i] );	 	
	 }	

	 for(int i=0; i<lftx.size(); i++){
		 	 w[lfty[i]][lftx[i]] = func_wbnd(1, 0, lftx[i], lfty[i] );	 	
	 }	
}

void Fluid::setbnd_psi(){}
bool vs = true;
void Fluid::initialize(float fx[M][N], float fy[M][N], float p[M][N], vector<int> Bx, vector<int> By, vector<float> vBx, vector<float> vBy ){

	if(vs){
		init_boundary(Bx, By, vBx, vBy);
		init_v(fx, fy);
	  init_w();
		setbnd_v();
		setbnd_w();
	}
	else{
		init_boundary(Bx, By, vBx, vBy);	
		init_v(fx, fy);
		setbnd_v();
	}
	
	cout<<"Parameters in dimentionless units"<<endl;	
	cout<<"dx = "<<DX<<", dy = "<<DY<<", dt = "<<dt<<endl;
  cout<<X1<<" <= x <= "<<X2;
  cout<<"  , "<<Y1<<" <= y <= "<<Y2<<endl;	
  cout<<"###############################"<<endl;
  cout<<"Reynolds number of flow : "<<Re<<endl;
	cout<<"Euler number of flow : "<<Eu<<endl;
  cout<<"###############################"<<endl;

}


void Fluid::vorticity_stream(){
	float vx0[M][N], vy0[M][N];
	
	poisson_solve(psi, w, psibx, psiby, psib );
	setbnd_w();	

//	predictor_corrector(w, vx, vy, sw, dt);
//	Explicit_FTCS_rk2(w , vx, vy, sw, dt);
	explicit_implicit(w, vx, vy, sw, bx, by, wb, dt );

	setbnd_w();	

	poisson_solve(psi, w, psibx, psiby, psib );

  ddx_cd( psi, vy0 );
  ddy_cd( psi, vx0 );

  for(int j=0; j<M; j++){
    for(int i=0; i<N; i++){
			vx[j][i] =  vx0[j][i]/(2.0*DY);
			vy[j][i] = -vy0[j][i]/(2.0*DX);
		}
	}
	setbnd_v();	
  cout<<Max(vx,vy)<<"	"<<t<<"		"<<Max(w)<<endl;	
	t = t+dt;

}


void Fluid::pressure_corrector(){
	
	float dvx[M][N], dvy[M][N], dxp[M][N]={}, dyp[M][N]={};
	float pc[M][N], dxpc[M][N], dypc[M][N];
	float cf = 2*DX*DX*DY*DY/( (DX*DX+DY*DY)*(2.0*dt*Eu) ), div, error;
	float vx0[M][N], vy0[M][N];
	
	vector<int> tx, ty;
	vector<float> p0;

	for(int i=0; i<bx.size(); i++){
		p0.push_back(0);
	}

	vort();
	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			vx0[j][i] = vx[j][i];
			vy0[j][i] = vy[j][i]; 
		}
	}	
	
	for(int itr=0; itr<10; itr++){

		for(int j=0; j<M; j++){
			for(int i=0; i<N; i++){
				vx[j][i] = vx0[j][i];
				vy[j][i] = vy0[j][i]; 
			}
		}	

		ddx_cd( P, dxp );
		ddy_cd( P, dyp );
		
  	for(int j=0; j<M; j++){
  		for(int i=0; i<N; i++){
				dxp[j][i] = Eu*dxp[j][i]/(2.0*DX) - svx[j][i] ;
				dyp[j][i] = Eu*dyp[j][i]/(2.0*DY) - svy[j][i] ; 
			}
		}	

		//	Explicit_FTCS_rk2
		predictor_corrector(vx, vx0, vy0, dxp, bx, by, vxb, dt );
		predictor_corrector(vy, vx0, vy0, dyp, bx, by, vyb, dt );

		ddx_cd( vx, dvx );
	  ddy_cd( vy, dvy ); 

	  for(int j=0; j<M; j++){
		  for(int i=0; i<N; i++){
				pc[j][i] =  -cf*( dvx[j][i]/(2.0*DX)+ dvy[j][i]/(2.0*DY));
			}
		}
		
//	  for(int j=0; j<M; j++){
//		  for(int i=0; i<N; i++){
//				dvx[j][i] =  ( dvx[j][i]/(2.0*DX)+ dvy[j][i]/(2.0*DY))/(Eu*dt);
//			}
//		}
//		poisson_solve(pc, dvx, bx, by, p0 );


		ddx_cd( pc, dxpc );
		ddy_cd( pc, dypc );
//		cout<<error<<"	"<<Max(vx,vy)<<"	"<<Max(dxpc, dypc)<<"	"<<Max(P)<<endl;
//		cout<<tx.size()<<endl;

	  for(int j=0; j<M; j++){
		  for(int i=0; i<N; i++){
				vx[j][i] = vx[j][i] - Eu*dt*dxpc[j][i]/(2.0*DX);
				vy[j][i] = vy[j][i] - Eu*dt*dypc[j][i]/(2.0*DY);
				P[j][i] = P[j][i] + pc[j][i];	
			}
		}

		ddx_cd( vx, dvx );
	  ddy_cd( vy, dvy );

	  error =0;
	  float max;	

		for(int j=M/2-M/4+2; j< M/2+M/4-2; j++){
			for(int i=N/2-N/4+2; i< N/2+N/4-2; i++){
  			div = abs(dvx[j][i]/(2.0*DX) + dvy[j][i]/(2.0*DY));
      	error = error + div; 
			}
		}
		setbnd_v();

		error = abs(error)/(M*N);	
  // 	cout<<itr<<" error : "<<error<<"		"<<max<<endl;    
    if(error < 0.01) {break;}

  } //pressure correction iteration
  //cout<<"Ffffffffffffffffffff "<<Max(vx,vy)<<"	"<<t<<"		"<<error<<endl;	  
	t = t + dt;

}

void Fluid::write_to_file(int x1, int y1, int x2, int y2, float t1, float t2, char filename[]){
  
  ofstream snap_t1, snap_t2, bnd;
  const int len = strlen(filename) + 25; 
	char st1[len], st2[len], bndr[len];
	int nop = (t2-t1)/dt+1;
	
	sprintf( st1, "%s_initial.txt", filename);
	sprintf( st2, "%s_final.txt", filename);
  sprintf( bndr, "%s_boundary.txt", filename);

	bnd.open(bndr); 
  for(int i=0; i<by.size(); i++ ){
    bnd<<"		"<<bx[i]<<"		"<<by[i]<<endl;    
  }  
  bnd.close();

  float t=t1;  
  for(long int frame=0; frame<nop; frame++){
		cout<<"t=	"<<t<<endl;
  	if(frame==0){
			snap_t1.open(st1);
		  for(int j=y1; j<=y2; j++){
		    for(int i=x1; i<=x2; i++){
	snap_t1<<"    "<<frame<<"		"<<t<<"   "<<i<<"   "<<j<<"   "<<vx[j][i]<<"   "<<vy[j][i]<<"   "<<P[j][i]<<"	"<<w[j][i]<<"	"<<psi[j][i]<<endl; 
		    } //x
		  }  //y
		  snap_t1.close();
  	}

  	if(frame==nop-1){
  		snap_t2.open(st2);
		  for(int j=y1; j<=y2; j++){
		    for(int i=x1; i<=x2; i++){
	snap_t2<<"    "<<frame<<"		"<<t<<"   "<<i<<"   "<<j<<"   "<<vx[j][i]<<"   "<<vy[j][i]<<"   "<<P[j][i]<<"	"<<w[j][i]<<"	"<<psi[j][i]<<endl; 
		    } //x
		  }  //y
			snap_t2.close();
  	}
 //   Evolve();
//   	vorticity_stream();
 	pressure_corrector();
    t = t1 + (frame+1)*dt; 
  } //frame  
}

void read_from_file(float &T, float fx[M][N], float fy[M][N], float p[M][N], float w[M][N], float psi[M][N], vector<int> &Bx, vector<int> &By, char filename[]){

	float VX, VY, TP , BX, BY, W, PSI;
	int frame,x,y;
	const int len = strlen(filename) + 25;   
	char snap[len], bndr[len];
	sprintf( snap, "%s_final.txt", filename);
  sprintf( bndr, "%s_boundary.txt", filename);
  ifstream input, bnd;
	input.open(snap);
  for(int j=0; j<M; j++){
    for(int i=0; i<N; i++){
       input>>frame>>T>>x>>y>>VX>>VY>>TP>>W>>PSI;
       fx[j][i]=VX;
       fy[j][i]=VY;
       p[j][i]=TP; 
       w[j][i]=W; 
       psi[j][i]=PSI; 
    } //x
  }  //y
	input.close();
	
	bnd.open(bndr);
	while(!bnd.eof()){
		bnd>>BX>>BY;
		Bx.push_back(BX);
		By.push_back(BY);	
	}
}













