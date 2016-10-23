//#include <mkl_lapacke.h> 
#include <lapacke.h> 
/*If you want to use lapacke then use   #include <lapacke.h> */
class matrix
{
	public:
	void diag(double *  a,double *  w,int n ) ; 
	void copy(double *  a,double *  b,int n ) ; 
	void dispMatrix(double *  a,int n ) ; 
	void checkHar(double *  a,int n ) ; 
	void dispArray(double *  a,int n ) ; 
	void ArrayCopy(double *  a,double *  b,int n ) ; 
	void copyDiag(double *  a,double *  b) ; 
} matrix ;
void matrix::ArrayCopy(double *  a,double *  b,int n)
{
	for(long int i=0;i<n;i++)
		a[i]=b[i] ; 
}
void matrix::dispArray(double *  a,int n)
{
	
	//for(long int i=0;i<n;i++)
	{
	for(long int j=0;j<n;j++)
		std::cout<<a[j]<<"\t"  ; 
	std::cout<<std::endl;
	}
}	
void matrix::checkHar(double *  a,int n)
{
	
	for(long int i=0;i<n;i++)
	for(long int j=i+1;j<n;j++)
	if (abs(a[i*n+j]-a[i+n*j]) > 1e-5 )
		std::cout<<"Error : Not Harmitian"<<std::endl;
}	
void matrix::dispMatrix(double *  a,int n)
{
	
	for(long int i=48;i<64;i++)
	{
	for(long int j=16;j<32;j++)
		std::cout<<int(100*a[i*n+j])<<"\t"  ; 
	std::cout<<std::endl;
	}
}	
void matrix::copy(double *  a,double *  b,int n)
{
	for(long int i=0;i<n*n;i++)
		a[i]=b[i] ; 
}
void matrix::diag(double *  a,double *  w,int n)
{
	int info= LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, n, w );
        if( info > 0 ) 
		std::cout<<"Not able to solve Eigen value problem."<<std::endl;
}
