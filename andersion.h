#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <algorithm>  
class QDanderson
{
    public:
	void initMatrix( double *  a ) ;// 
	void testDelete( double *  a ) ;// 
	void  addaSite(double * a,double * b, double * c ); 
	void  calcSpec(double * a,double * b,std::ofstream& pfile ) ;// 
}QDanderson ;
class param
{
    public:
	double U,B,eps,V,LAMBDA,sqLAMBDA ; 
	int startDim,wchain,bathDim,currDim,maxStates,dotDim,noStates,preDim;
	void setparam() ;
	void rescale() ;
	double hopping(int site) ;
}param ;
void param::setparam( )
{
	maxStates=5000;
	dotDim=4;
	wchain=1;
	bathDim=std::pow(4,wchain);
	startDim=dotDim*bathDim;
	currDim=dotDim*bathDim;
	noStates=dotDim*bathDim;
	preDim=dotDim*bathDim;
	
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
	return (1.0+1.0/LAMBDA)*(1-std::pow(LAMBDA,-site-1))*0.50/std::sqrt((1.0-std::pow(LAMBDA,-2.*site-1))*(1-std::pow(LAMBDA,-2.*site-3))) ;
}
void  QDanderson::calcSpec(double * a,double * b,std::ofstream& pfile ) 
{
	double tsup1,tsdw1,tsup2,tsdw2 ;
	int n=param.currDim;
	int m=param.currDim/4;
	for(int iw=1;iw<n;iw++){
		tsup1=0;tsdw1=0;tsup2=0;tsdw2=0;
		for( int nci=0;nci<m; nci++){
		//Up spin
		tsup1+= (a[iw*n+(1+4*nci)]*a[0*n+(0+4*nci)] + a[iw*n+(3+4*nci)]*a[0*n+(2+4*nci)] ) ; 
		tsdw1+= (a[iw*n+(2+4*nci)]*a[0*n+(0+4*nci)] - a[iw*n+(3+4*nci)]*a[0*n+(1+4*nci)] ) ; 
		tsup2+= (a[0*n+(1+4*nci)]*a[iw*n+(0+4*nci)] + a[0*n+(3+4*nci)]*a[iw*n+(2+4*nci)] ) ; 
		tsdw2+= (a[0*n+(2+4*nci)]*a[iw*n+(0+4*nci)] - a[0*n+(3+4*nci)]*a[iw*n+(1+4*nci)] ) ; 
		}
	std::cout<<"----:  "<<b[iw]<<std::endl;
	pfile<<(b[iw]-b[0])<<" "<<std::pow(tsup1,2)<< " "<<std::pow(tsdw1,2)<<std::endl;
	pfile<<(b[0]-b[iw])<<" "<<std::pow(tsup2,2)<< " "<<std::pow(tsdw2,2)<<std::endl;
	}
}
void  QDanderson::addaSite(double * a,double * b, double * c ) 
{
	double aa = param.hopping(param.wchain) ;
	int m=param.currDim;
	int n=param.currDim/4;
	int p=param.noStates;
	int prem=param.preDim; 
	int rs=param.preDim/4; 

		
		
	std::cout<<"matrix dim m: "<<m<<"n: "<<n<<"p: "<<p <<std::endl;
	std::cout<<"Previous Matrix: "<<prem<<std::endl;
	std::cout<<"Hopping element : " <<aa <<std::endl;
	

	for(long int i=0;i<p;i++){
		a[i*m+i]=param.sqLAMBDA*c[i] ; 
		a[(i+p)*m+(i+p)]=param.sqLAMBDA*c[i] ; 
		a[(i+2*p)*m+(i+2*p)]=param.sqLAMBDA*c[i] ; 
		a[(i+3*p)*m+(i+3*p)]=param.sqLAMBDA*c[i] ; 
	}

	for(int i=0;i<p;i++)
	for(int j=0;j<p;j++)
	{
		double	upn=0,dwn=0;
		for( int nci=0;nci<prem/4; nci++){
		upn += aa*(b[i*prem+(nci+rs*1)]*b[j*prem+(nci+rs*0)] + b[i*prem+(nci+rs*3)]*b[j*prem+(nci+rs*2)]) ; 
		dwn += aa*(b[i*prem+(nci+rs*2)]*b[j*prem+(nci+rs*0)] - b[i*prem+(nci+rs*3)]*b[j*prem+(nci+rs*1)]) ;
		}
		//Up spin
		a[(i+n*0)*m+(j+n*1)] += upn; 
		a[(i+n*2)*m+(j+n*3)] -= upn;
		a[(i+n*0)+m*(j+n*1)] += upn; 
		a[(i+n*2)+m*(j+n*3)] -= upn;

		//Down spin
		a[(i+n*0)*m+(j+n*2)] += dwn;
		a[(i+n*1)*m+(j+n*3)] += dwn;
		a[(i+n*0)+m*(j+n*2)] += dwn;
		a[(i+n*1)+m*(j+n*3)] += dwn;
	}

}	
void QDanderson::initMatrix(double * a)
{

	//Dot Hamiltonian H = eps*n_d - B/2 (n_{d,up}-n_{d,down}) + U . n_{d,up} . n_{d,down}
	//Also Rescale the Dot hamitonian
	int dotDim=4;
	double rescale=std::pow(param.LAMBDA,-0.5) ;
	std::cout<<"Rescale : "<<rescale<<std::endl	;
	for(int nci=0;nci<param.bathDim; nci++)
	{
	a[(1+dotDim*nci)*param.startDim+(1+dotDim*nci)] = (param.eps - 1./2. * param.B)*rescale ; 
	a[(2+dotDim*nci)*param.startDim+(2+dotDim*nci)] = (param.eps + 1./2. * param.B)*rescale ; 
	a[(3+dotDim*nci)*param.startDim+(3+dotDim*nci)] = (2.*param.eps + param.U)*rescale ; 
	}
	//Connect the dot with the first lattice of the wilson Chain:
	//d^\dag_{\up} c_{up}
	a[(1+dotDim*0)*param.startDim+(0+dotDim*1)]=param.V*rescale;
	a[(3+dotDim*0)*param.startDim+(2+dotDim*1)]=param.V*rescale;
	a[(1+dotDim*2)*param.startDim+(0+dotDim*3)]=-param.V*rescale;
	a[(3+dotDim*2)*param.startDim+(2+dotDim*3)]=-param.V*rescale;

	//d^\dag_{\down} c_{down}	
	a[(2+dotDim*0)*param.startDim+(0+dotDim*2)]=param.V*rescale;
	a[(3+dotDim*0)*param.startDim+(1+dotDim*2)]=-param.V*rescale;
	a[(2+dotDim*1)*param.startDim+(0+dotDim*3)]=param.V*rescale;
	a[(3+dotDim*1)*param.startDim+(1+dotDim*3)]=-param.V*rescale;
	
	//	Complex comjugate:
	//d^\dag_{\up} c_{up}
	a[(1+dotDim*0)+param.startDim*(0+dotDim*1)]=param.V*rescale;
	a[(3+dotDim*0)+param.startDim*(2+dotDim*1)]=param.V*rescale;
	a[(1+dotDim*2)+param.startDim*(0+dotDim*3)]=-param.V*rescale;
	a[(3+dotDim*2)+param.startDim*(2+dotDim*3)]=-param.V*rescale;

	//d^\dag_{\down} c_{down}	
	a[(2+dotDim*0)+param.startDim*(0+dotDim*2)]=param.V*rescale;
	a[(3+dotDim*0)+param.startDim*(1+dotDim*2)]=-param.V*rescale;
	a[(2+dotDim*1)+param.startDim*(0+dotDim*3)]=param.V*rescale;
	a[(3+dotDim*1)+param.startDim*(1+dotDim*3)]=-param.V*rescale;
 
}

