#include <limits>
#include <cmath>
#include "Includes/MatrixInverse_Test.h"


bool inv_test(const std::vector< std::vector< double> >& matrix)
{
    std::vector< std::vector<double> > minv=inv(matrix);

    std::vector< std::vector<double> > res=matrix;

    double max=0;
    std::size_t n=matrix.size();
    for (std::size_t i=0; i<n; ++i)
	for (std::size_t j=0; j<n; ++j)
	    if (std::abs(matrix[i][j])>max)
		max=std::abs(matrix[i][j]);

    bool result=true;
    double tolerance=sqrt(std::numeric_limits<double>::epsilon())*max;

    for (std::size_t i=0; i<n; ++i)
	for (std::size_t j=0; j<n; ++j)
	{
	    res[i][j]=0;
	    for (std::size_t k=0; k<n; ++k)
		res[i][j]+=matrix[i][k]*minv[k][j];
	    if (i!=j)
	    {if (std::abs(res[i][j])>tolerance)
		    result=false;
	    }
	    else
	    {
		    if (std::abs(res[i][j]-1.0)>tolerance)
			result=false;
	    }
	}



    std::cout<<"input matrix"<<matrix;
    std::cout<<"\n----------------------------------------------------------\n";
    std::cout<<"inverse"<<minv;
    std::cout<<"\n----------------------------------------------------------\n";

    if (result)
    std::cout<<"\n TEST PASSED !!!\n";
    else
	std::cout<<"\n TEST FAIL !!!\n";


    std::cout<<"matrix * inverse="<<res;

    return result;


}


bool inv_test()
{
    std::vector< std::vector <double> > A(3,std::vector<double>(3));
    A[0][0]=12;
    A[0][1]=1;
    A[0][2]=-14;
    A[1][0]=2;
    A[1][1]=41;
    A[1][2]=99;
    A[2][0]=-0.02;
    A[2][1]=410;
    A[2][2]=-999;

    return inv_test(A);

}



bool det_test()
{
    std::vector< std::vector <double> > A(3,std::vector<double>(3));
    A[0][0]=14;
    A[0][1]=37;
    A[0][2]=44;
    A[1][0]=37;
    A[1][1]=110;
    A[1][2]=127;
    A[2][0]=44;
    A[2][1]=127;
    A[2][2]=149;

    double d=det(A);
    std::cout<<"el det es "<<d;
    return d==225;

}
