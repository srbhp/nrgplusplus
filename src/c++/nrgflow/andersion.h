#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <algorithm>  
#include <mkl.h>
#include <omp.h>
#include "matrix.h"
class QDanderson
{
    public:
	void initMatrix( double *  a ) ;// 
	double spheat( double *  w ) ;// 
	void shiftEnergy( double *  w ,int n) ;// 
	void addaSite(double * a,double * b, double * c ,double rescale); 
	void calcSpec(double * a,double * b,std::ofstream& pfile ) ;// 
	void updateUP(double *a ,double * up,double * oldup);
	void initOperator(double *a, double * up);

}QDanderson ;
class param
{
    public:
	double U,B,eps,V,LAMBDA,sqLAMBDA ; 
	int wchain,currDim,maxStates,noStates,preDim;
	void setparam() ;
	void rescale() ;
	double hopping(int site) ;
}param ;
void param::setparam( )
{
	maxStates=5000;
	wchain=-1;
	currDim=4;
	noStates=4;
	preDim=4;
	
	U=0.10 ;
	V=0.250;  
	B=0.0 ;
	LAMBDA=2.50;
	sqLAMBDA=std::sqrt(LAMBDA);
	eps = -U/2.0 ;
}
void param::rescale( )
{
	U=U*2.0/(1.0 +1.0/LAMBDA ) ;
	eps=eps*2.0/(1.0 +1.0/LAMBDA ) ;
	B=B*2.0/(1.0 +1.0/LAMBDA ) ;
}

double  param::hopping(int site )
{
	return  (1.0+1.0/LAMBDA)*(1-std::pow(LAMBDA,-site-1))*0.50/std::sqrt((1.0-std::pow(LAMBDA,-2.*site-1))*(1-std::pow(LAMBDA,-2.*site-3))) ;
}
double QDanderson::spheat( double *  w)
{
	double a1=0.0 ; 
	for(int i=0;i<param.currDim;i++)
	a1=a1+ w[i]*w[i]*std::exp(-w[i]) ;
	return a1 ; 
}
void  QDanderson::shiftEnergy( double *  w ,int n)
{
	double * bs=new double[n]() ; 
	for(long int i=0;i<n;i++)
		bs[i]=w[i];
	for(long int i=0;i<n;i++)
		w[i]=bs[i] - bs[0] ;

	delete [] bs ; 
}
void  QDanderson::calcSpec(double * a,double * b,std::ofstream& pfile ) 
{
	int n=param.currDim;
	int m=param.currDim/4;
	double a1=0,a2=0;

	/*
	for(int iw=1;iw<n;iw++){
	a1=a1+std::pow(a[0+n*iw],2);
	a2=a2+std::pow(a[0*n+iw],2);
	}*/
	
	for(int iw=0;iw<n;iw++){
	a1=a1+a[iw*n+iw] ;
	a2=a2+std::exp(-b[iw]) ;
	}
	pfile<<2.0*std::pow(param.LAMBDA,-param.wchain/2.0)<<" "<<a1<<" "<<a2 <<std::endl;
}
void  QDanderson::addaSite(double * a,double * b, double * c,double aa ) 
{
	double diaa =param.sqLAMBDA ;
	if( param.wchain < -0.5 ){
		// Parameters for the first Wilson's site 
		aa=param.V*std::pow(param.LAMBDA,-0.5)  ; 
		diaa =1.0 ;
	}

	int m=param.currDim;
	int n=param.currDim/4;
	int p=param.noStates;
	int prem=param.preDim; 
	int rs=param.preDim/4; 

	std::cout<<"Wilson's chain: "<<param.wchain<<" hopping int. "<<aa<<" Rescale: "<<diaa <<std::endl;
	std::cout<<"matrix dim m: "<<m<<" n: "<<n<<" p: "<<p <<std::endl;
	std::cout<<"Previous Matrix: "<<prem<<std::endl;
	
	#pragma omp parallel for
	for(long int i=0;i<p;i++){
		a[i*m+i]=diaa*c[i] ; 
		a[(i+p)*m+(i+p)]=diaa*c[i] ; 
		a[(i+2*p)*m+(i+2*p)]=diaa*c[i] ; 
		a[(i+3*p)*m+(i+3*p)]=diaa*c[i] ; 
	}

	#pragma omp parallel for  collapse(2)
	for(int i=0;i<p;i++)
	for(int j=0;j<p;j++)
	{
		double	upn=0,dwn=0;
		double	upd=0,dwd=0;
		for(int nci=0;nci<prem/4; nci++){
		upn += (b[i*prem+(nci+rs*1)]*b[j*prem+(nci+rs*0)] + b[i*prem+(nci+rs*3)]*b[j*prem+(nci+rs*2)]) ; 
		dwn += (b[i*prem+(nci+rs*2)]*b[j*prem+(nci+rs*0)] - b[i*prem+(nci+rs*3)]*b[j*prem+(nci+rs*1)]) ;

		upd += (b[i*prem+(nci+rs*0)]*b[j*prem+(nci+rs*1)] + b[i*prem+(nci+rs*2)]*b[j*prem+(nci+rs*3)]) ; 
		dwd += (b[i*prem+(nci+rs*0)]*b[j*prem+(nci+rs*2)] - b[i*prem+(nci+rs*1)]*b[j*prem+(nci+rs*3)]) ;

		/*
		upn += (b[i+prem*(nci+rs*1)]*b[j+prem*(nci+rs*0)] + b[i+prem*(nci+rs*3)]*b[j+prem*(nci+rs*2)]) ; 
		dwn += (b[i+prem*(nci+rs*2)]*b[j+prem*(nci+rs*0)] - b[i+prem*(nci+rs*3)]*b[j+prem*(nci+rs*1)]) ;

		upd += (b[i+prem*(nci+rs*0)]*b[j+prem*(nci+rs*1)] + b[i+prem*(nci+rs*2)]*b[j+prem*(nci+rs*3)]) ; 
		dwd += (b[i+prem*(nci+rs*0)]*b[j+prem*(nci+rs*2)] - b[i+prem*(nci+rs*1)]*b[j+prem*(nci+rs*3)]) ;
		*/
		}
		//Up spin
		a[(i+p*0)*m+(j+p*1)] = upn*aa; 
		a[(i+p*2)*m+(j+p*3)] = -upn*aa;
		a[(i+p*1)*m+(j+p*0)] = upd*aa; 
		a[(i+p*3)*m+(j+p*2)] = -upd*aa;

		//Down spin
		a[(i+p*0)*m+(j+p*2)] = dwn*aa;
		a[(i+p*1)*m+(j+p*3)] = dwn*aa;
		a[(i+p*2)*m+(j+p*0)] = dwd*aa;
		a[(i+p*3)*m+(j+p*1)] = dwd*aa;
	}
	
	
}	
void QDanderson::initMatrix(double * a)
{
	std::cout<<"eps "<<param.eps<<" U "<<param.U<<std::endl;
	//Dot Hamiltonian H = eps*n_d - B/2 (n_{d,up}-n_{d,down}) + U . n_{d,up} . n_{d,down}
	//Also Rescale the Dot hamitonian
	int dotDim=4;
	double rescale=1.0/std::pow(param.LAMBDA,0.5) ;
	std::cout<<"Rescale : "<<rescale<<std::endl	;
	a[(1)*dotDim+(1)] = (param.eps - 1./2. * param.B)*rescale ; 
	a[(2)*dotDim+(2)] = (param.eps + 1./2. * param.B)*rescale ; 
	a[(3)*dotDim+(3)] = (2.*param.eps + param.U)*rescale ; 
	//Connect the dot with the first lattice of the wilson Chain:
	//d^\dag_{\up} c_{up}

}

void QDanderson::initOperator(double *a, double * up)
{
	int m=param.currDim;
	for(long int i=0;i<m;i++)
	for(long int j=0;j<m;j++)
	{
		//up[i*m+j] = a[i+m*1]*a[j+m*0] + a[i+m*3]*a[j+m*2 ] ;
		up[i*m+j] = a[i+m*3]*a[j+m*3]  ;
	}

}

void QDanderson::updateUP(double *a ,double * up,double * oldup)
{
	int long m=param.currDim;
	int long p=m/4;
	std::cout<<"Updating the Operators .. " << std::endl;
	std::cout<<"Old: "<<p <<" "<<m<<std::endl;
	double * uptm = new double [m*m] ();
	double * tmat = new double [m*m] ();
	//if( !uptm == nullptr || tmat == nullptr )
	if( !uptm  || !tmat )
		std::cout<<"Unable to allocate memory "<<std::endl;
		
	#pragma omp parallel for  collapse(2)
	for(long int i=0;i<p;i++)
	for(long int j=0;j<p;j++){
		uptm[i*m+j]=oldup[i*p+j] ; 
		uptm[(i+p)*m+(j+p)]=oldup[i*p+j]  ; 
		uptm[(i+2*p)*m+(j+2*p)]=oldup[i*p+j] ; 
		uptm[(i+3*p)*m+(j+3*p)]=oldup[i*p+j] ; 
	}
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
                m, m, m, 1, a, m, uptm, m, 0, tmat, m);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
                m, m, m, 1, tmat, m, a, m, 0, up, m);
	
	//matrix.dispMatrix(uptm,m) ;
	delete [] uptm;
	delete [] tmat;
	
	
	/*
	#pragma omp parallel for  collapse(2)
	for(long int in=0;in<m;in++)
	for(long int jn=0;jn<m;jn++)
	{	
		for(long int i=0;i<p;i++)
		for(long int j=0;j<p;j++)
		for(long int k=0;k<4;k++)
		{
		
			up[in*m+jn] += a[in+m*(i*4+k )]*oldup[i*p+j]*a[jn+m*(j*4+k)] ;
		}
	}
	*/
}

