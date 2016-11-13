#include <mkl_lapacke.h> 
//#include <lapacke.h> 
/*If you want to use lapacke then use   #include <lapacke.h> */
class matrix
{
	public:
	void diag(double *  a,double *  w,int n ) ; 
	void copyMatrix(double *  a,double *  b,int n ) ; 
	void dispMatrix(double *  a,int n ) ; 
	void dispDiag(double *  a,int n ) ; 
	void checkHar(double *  a,int n ) ; 
	void dispArray(double *  a,int n ) ; 
	void copyArrG(double *  a,double *  b,int n ) ; 
} matrix ;
void matrix::copyArrG(double *  a,double *  b,int n)
{
	double ats=b[0] ;
	#pragma omp parallel for
	for(long int i=0;i<n;i++)
		a[i]=b[i]-ats; 
}
void matrix::dispArray(double *  a,int n)
{
	
	for(long int j=0;j<n;j++)
		std::cout<<a[j]<<"\t"  ; 
	std::cout<<std::endl;
}	
void matrix::checkHar(double *  a,int n)
{
	
	for(long int i=0;i<n;i++)
	for(long int j=i+1;j<n;j++)
	if (abs(a[i*n+j]-a[i+n*j]) > 1e-5 )
		std::cout<<"Error : Not Harmitian"<<std::endl;
}	
void matrix::dispDiag(double *  a,int n)
{
	
	//for(long int i=0;i<n;i++)
	{
	for(long int j=0;j<n;j++)
		std::cout<<(a[j*n+j])<<"\t"  ; 
	std::cout<<std::endl;
	}
}	
void matrix::dispMatrix(double *  a,int n)
{
	
	for(long int i=0;i<n;i++)
	{
	for(long int j=0;j<n;j++)
		std::cout<<a[i*n+j]<<"\t"  ; 
	std::cout<<std::endl;
	}
}	
void matrix::copyMatrix(double *  a,double *  b,int n)
{
	#pragma omp parallel for
	for(long int i=0;i<n*n;i++)
		a[i]=b[i] ; 
}
void matrix::diag(double *  a,double *  w,int n)
{
	int info= LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, n, w );
        if( info > 0 ) 
		std::cout<<"Not able to solve Eigen value problem."<<std::endl;
}
