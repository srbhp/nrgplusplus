#include "andersion.h"
#include "matrix.h"
//#include <chkp.h>
int main()
{
	std::cout<<"Starting NRG .... "<<std::endl;
	param.setparam() ;
	param.maxStates=500;
	param.U=1e-3 ; 
	param.V=0.004 ; 
	param.eps=-param.U/2.0 ;
	std::ostringstream ss;
	ss <<"spectr-U"<<param.U<<".dat";
	std::string s(ss.str());
	std::ofstream pfile(s) ; 

	//Create first Initial Matrix 
	double * a = new double [param.startDim*param.startDim] () ; //matrix
	double * exa = new double [param.startDim*param.startDim] () ; //matrix
	double * w = new double [param.startDim] () ; //eigen value
	double * exw = new double [param.startDim] () ; //eigen value
	QDanderson.initMatrix(a) ; //create matrix
	std::cout<<"Starting NRG Dim: "<<param.startDim<<std::endl;
	

	matrix.checkHar(a,param.currDim) ; 
	matrix.diag(a,w,param.startDim) ; //Diagonalize the matrix
	matrix.dispArray(w,param.startDim) ;

	int Noiter=5,iter=1;
	for(int iter=1; iter<Noiter;iter++)
	{
	

	std::cout<<"Starting Wilson's chain: "<<param.wchain<<std::endl;
	std::cout<<"Starting Total Dim: "<<param.currDim<<std::endl;
	//Add a site to the wilson Chain.
	param.currDim=param.noStates*4 ; 
	delete [] exa ; 
	exa = new (std::nothrow) double [param.currDim*param.currDim]() ; 
	if (!exa) std::cout<<"Unable to allocate memory for exa"<<std::endl;  
	QDanderson.addaSite(exa,a,w) ;
	
	
	delete [] exw; 
	exw = new (std::nothrow) double [param.currDim] () ; //eigen value
	if (!exw) std::cout<<"Unable to allocate memory for exw"<<std::endl;  
	
	matrix.diag(exa,exw,param.currDim) ;
	for(int i =0;i<10;i++)
		std::cout<<exw[i]<<" "<<w[i]<<std::endl; 
	for(int i =param.currDim-10;i<param.currDim;i++)
		std::cout<<exw[i]<<" "<<w[i]<<std::endl ;
	
	delete [] a;
	a = new (std::nothrow)  double [param.currDim*param.currDim] () ; //matrix
	if (!a) std::cout<<"Unable to allocate memory for a"<<std::endl;  
	delete [] w ;
	w = new (std::nothrow)  double [param.currDim] () ; //matrix
	if (!w) std::cout<<"Unable to allocate memory for w"<<std::endl;  
	matrix.copy(a,exa,param.currDim) ; //Keep a copy of the matrix 
	matrix.ArrayCopy(w,exw,param.currDim) ; //Keep a copy of the matrix 
	
	std::cout<<" -----------------------------------------------"<<std::endl;
	for(int i =0;i<10;i++)
		std::cout<<exw[i]<<" "<<w[i]<<std::endl; 
	for(int i =param.currDim-10;i<param.currDim;i++)
		std::cout<<exw[i]<<" "<<w[i]<<std::endl ;
	
	param.preDim=param.currDim	;
	param.wchain=param.wchain+1;
	if(4*param.noStates < param.maxStates )
		param.noStates = 4*param.noStates ; 
	else 
		param.noStates = param.maxStates ;


	}
	std::cout<<" -------------iFINAL----------------------------------"<<std::endl;
	for(int i =0;i<10;i++)
		std::cout<<exw[i]<<" "<<w[i]<<std::endl; 
	for(int i =param.currDim-10;i<param.currDim;i++)
		std::cout<<exw[i]<<" "<<w[i]<<std::endl ;
	
	std::cout<<"Calc for spectrum started .. "<<std::endl;
	QDanderson.calcSpec(exa,exw,pfile ) ;

	
	
	pfile.close();
	delete [] a;
	delete [] w;
	delete [] exa;
	delete [] exw;		
	return 0;

}

