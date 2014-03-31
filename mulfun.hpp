#include <vector>
#include <linbox/matrix/blas-matrix.h>
#include "linbox/algorithms/blas-domain.h"


template<class Field>
void karatsuba(const Field& F,
		std::vector<LinBox::BlasMatrix<Field> >& C,
		const std::vector<LinBox::BlasMatrix<Field> >& A,
		const std::vector<LinBox::BlasMatrix<Field> >& B);


template<class Field>
void efmul(const Field& F,
		std::vector<LinBox::BlasMatrix<Field> >& C,
		const std::vector<LinBox::BlasMatrix<Field> >& A,
		const std::vector<LinBox::BlasMatrix<Field> >& B,
		const std::vector<typename Field::Element>& P)
{
	
	std::vector<LinBox::BlasMatrix<Field> > tempC (A.size()+B.size()-1,
			LinBox::BlasMatrix<Field> (F, A[0].rowdim(), B[0].coldim()));

	C = std::vector<LinBox::BlasMatrix<Field> >(P.size()-1,
				LinBox::BlasMatrix<Field>(F,A[0].rowdim(),B[0].coldim()));

	LinBox::BlasMatrixDomain<Field> BMD(F);
	LinBox::MatrixDomain<Field> MD(F);
	karatsuba(F,tempC,A,B);
	

	
	unsigned int i=tempC.size()-1;
	while(i >= P.size()-1)
	{
		unsigned int j=i+1-P.size();
		typename Field::Element c;
		LinBox::BlasMatrix<Field> tempM(F,A[0].rowdim(),B[0].coldim());
		for(unsigned int k=0; k<P.size(); ++k)
		{
			F.neg(c,P[k]);
			MD.mul(tempM,tempC[i],c);
			BMD.addin(tempC[j+k],tempM);
		}
		--i;
	}

	for(unsigned int i=0; i < C.size(); ++i)
	{
		C[i]=tempC[i];
	}
}

template<class Field>
void karatsuba(const Field& F,
		std::vector<LinBox::BlasMatrix<Field> >& C,
		const std::vector<LinBox::BlasMatrix<Field> >& A,
		const std::vector<LinBox::BlasMatrix<Field> >& B)
{
	const unsigned int n = A.size()/2;
	C = std::vector<LinBox::BlasMatrix<Field> >(A.size()+B.size()-1 , LinBox::BlasMatrix<Field>(F, A[0].rowdim(),B[0].coldim()));
	LinBox::BlasMatrixDomain<Field> BMD(F);
	if(n == 0) {
		BMD.mul(C[0],A[0],B[0]);
	}
	else {

		std::vector<LinBox::BlasMatrix<Field> > A1(A.begin(), A.begin() + n);
		std::vector<LinBox::BlasMatrix<Field> > A2(A.begin() + n, A.end());
		std::vector<LinBox::BlasMatrix<Field> > B1(B.begin(), B.begin() + n);
		std::vector<LinBox::BlasMatrix<Field> > B2(B.begin() + n, B.end());

		std::vector<LinBox::BlasMatrix<Field> > C1, C2, C3;
		karatsuba(F,C1,A1,B1);

		karatsuba(F,C2,A2,B2);

		for(unsigned int i = 0; i < n; ++i) {
		  BMD.addin(A1[i],A2[i]);
		  BMD.addin(B1[i],B2[i]);
		}
		karatsuba(F,C3,A1,B1);
		for(unsigned int i = 0; i < C1.size(); ++i) {
		  BMD.addin(C[i],C1[i]);
		  BMD.subin(C[i+n],C1[i]);
		  BMD.addin(C[i+2*n],C2[i]);
		  BMD.subin(C[i+n],C2[i]);
		  BMD.addin(C[i+n],C3[i]);
		}

	}	
}


