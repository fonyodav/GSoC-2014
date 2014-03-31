#include <iostream>
#include <vector>
#include <linbox/field/modular.h>
#include <linbox/vector/stream.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/matrix/blas-matrix.h>
#include <fstream>
#include "mulfun.hpp"



int main()
{
	const int P=7;  			// characteristic of the field
	const int K=4;				// power of p in the order of the field
	const int M=2;				// row-dimension of A
	const int N=3;				// coloumn-dimension of A and row-dimension of B
	const int L=4;				// coloumn-dimension of B

	typedef LinBox::Modular<int> Field;
	typedef std::vector<Field::Element> Polynomial;

	Field F(P);
	Polynomial pol(K+1);
	LinBox::RandomDenseStream<Field, Polynomial> factory1 (F, K+1);
	factory1 >> pol;
	F.init(pol[K],1);

	std::vector<LinBox::BlasMatrix<Field> > A(K,LinBox::BlasMatrix<Field>(F,M,N));
	std::ifstream inputA("a.txt");
	for(int i=0; i<K; i++) {
		A[i].read(inputA);
		/*A[i].write(std::cout);
		std::cout << std::endl;*/
	}

	std::vector<LinBox::BlasMatrix<Field> > B(K,LinBox::BlasMatrix<Field>(F,N,L));
		std::ifstream inputB("b.txt");
		for(int i=0; i<K; i++) {
			B[i].read(inputB);
			/*B[i].write(std::cout);
			std::cout << std::endl;*/
		}



	std::cout<<"A:\n";
	for(int i=0; i<M; ++i) {
		for(int j=0; j<N; ++j) {
			for(int h=K-1;h>=0;--h) {
				F.write(std::cout,A[h][i][j]);
				if(h>1) {
					std::cout<<"*x^"<<h<<"+";
				}
				else if(h == 1) {
					std::cout<<"*x+";
				}
			}
			std::cout<<"\t";
		}
		std::cout<<std::endl;
	}

	std::cout<<"B:\n";
	for(int i=0; i<N; ++i) {
		for(int j=0; j<L; ++j) {
			for(int h=K-1; h>=0; --h) {
				F.write(std::cout,B[h][i][j]);
				if(h>1) {
					std::cout<<"*x^"<<h<<"+";
				}
				else if(h == 1) {
					std::cout<<"*x+";
				}
			}
			std::cout<<"\t";
		}
		std::cout<<std::endl;
	}

	std::cout<<"polynom:\n";
	for(int h=K; h>=0; --h) {
		F.write(std::cout,pol[h]);
		if(h>1) {
			std::cout<<"*x^"<<h<<"+";
		}
		else if(h == 1) {
			std::cout<<"*x+";
		}
	}
	std::cout << std::endl;


	std::vector<LinBox::BlasMatrix<Field> > C(K,LinBox::BlasMatrix<Field>(F,M,L));
	efmul(F,C,A,B,pol);


	std::cout<<"C:\n";
	for(int i=0; i<M; ++i) {
		for(int j=0; j<L; ++j) {
			for(int h=K-1; h>=0; --h) {
				F.write(std::cout,C[h][i][j]);
				if(h>1) {
					std::cout<<"*x^"<<h<<"+";
				}
				else if(h == 1) {
					std::cout<<"*x+";
				}
			}
			std::cout<<"\t";
		}
		std::cout<<std::endl;
	}


}
