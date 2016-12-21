#include <mkl_lapacke.h> 
//#include <lapacke.h>

/*If you want to use lapacke then use   #include <lapacke.h> */
#include<iomanip>
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
	void chfermionOP(double * a, int n);
} matrix ;
void matrix::chfermionOP(double * a, int n)
{
	double a1=0,a2=0;
	double ct1=omp_get_wtime();

	int long m=n;


	double * tmat2 = new double [m*m] ();
	double * tmat = new double [m*m] ();
	//if( !uptm == nullptr || tmat == nullptr )
	if( !tmat2  || !tmat )
		std::cout<<"Unable to allocate memory "<<std::endl;

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
                m, m, m, 1, a, m, a, m, 0, tmat, m);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                m, m, m, 1, a, m, a, m, 0, tmat2, m);
	
	//matrix.dispMatrix(uptm,m) ;


	#pragma omp parallel for reduction(+:a1,a2)
	for(int i=0;i<n;i++){
		a1 = a1 + tmat[i*m+i] + tmat2[i*m+i] ; 
		for(int k=i+1;k<n;k++){
		a2=a2+tmat[i*m+k] + tmat2[i*m+k]+tmat[i+m*k] + tmat2[i+m*k] ;
	}
	}
	a1=a1/n ;
	a2=a2/n;
	std::cout<<"Sum of Dia and off-Dia elements are : "<<a1<<" "<<a2<<std::endl;
	a1=1.0/std::sqrt(a1) ; 

	#pragma omp parallel for
	for(long  int i=0;i<n*n;i++)
		a[i]=a[i]*a1 ;

	std::cout<<"Time spend to calc. anti-commu: "<<omp_get_wtime()-ct1<<"sec"<<std::endl;
	delete [] tmat;
	delete [] tmat2;
}

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
		std::cout<<std::fixed<<std::setprecision(3)<<a[i*n+j]<<" "  ; 
	std::cout<<std::endl;
	}
}	
void matrix::copyMatrix(double *  a,double *  b,int n)
{
	double ct1=	omp_get_wtime();
	#pragma omp parallel for
	for(long int i=0;i<n*n;i++)
		a[i]=b[i] ; 
	std::cout<<"Time spend to copy a Matrix: "<<omp_get_wtime()-ct1<<"sec"<<std::endl;
}
void matrix::diag(double *  a,double *  w,int n)
{
	double ct1=omp_get_wtime();
	int info= LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, n, w );
	if( info > 0 ) std::cout<<"Not able to solve Eigen value problem."<<std::endl;


	std::cout<<"Time spend to Diag: "<<omp_get_wtime()-ct1<<"sec"<<std::endl;

}
