#include <vector>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/matrix-domain.h>



template<class Field>
void efmul(const Field& F,
		std::vector<LinBox::BlasMatrix<Field> >& C,
		const std::vector<LinBox::BlasMatrix<Field> >& A,
		const std::vector<LinBox::BlasMatrix<Field> >& B,
		const std::vector<typename Field::Element>& P)
{
	std::vector<LinBox::BlasMatrix<Field> > tempC (A.size()+B.size()-1,
			LinBox::BlasMatrix<Field> (F, A[0].rowdim(), B[0].coldim()));

	std::cout << A[0].rowdim() << A[0].coldim() << B[0].rowdim() << B[0].coldim() << std::endl;

	C = std::vector<LinBox::BlasMatrix<Field> >(P.size()-1,
				LinBox::BlasMatrix<Field>(F,A[0].rowdim(),B[0].coldim()));

	LinBox::MatrixDomain<Field> MD(F);
	for(unsigned int i=0; i < A.size(); ++i) {
		for(unsigned int j=0; j < B.size(); ++j) {
			MD.axpyin(tempC[i+j],A[i],B[j]);
		}
	}

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
			MD.addin(tempC[j+k],tempM);
		}
		--i;
	}


	for(unsigned int i=0; i < C.size(); ++i)
	{
		C[i]=tempC[i];
	}

}
