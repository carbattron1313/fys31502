//
//  Project2.cpp
//
//  Elena Muñoz Rivas and Alejandro Carballido Mantecón.
//

#include "Project2.hpp"
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>

int main(){

    //To set up the tridiagonal A for N=6, we create the constats:
    int N = 6;
    int n = N+1;
    double h = 1./n;

    //Create de values of the vectors of the diagonal and triagonal for matrix A:
    double a = -1./(h*h);
    double d = 2./(h*h);
    
    //Create the matrices A, eigvec and the vector eigval:
    arma::mat A = create_tridiagonal(N,a,d,a);
    arma::vec eigval;
    arma::mat eigvec;

    //Use the armadillo eig_sym function:
    arma::eig_sym(eigval, eigvec, A);
 
    //Normalise the vector eigval:
    arma::vec eigval_n = arma::normalise(eigval,1);
    //Create the matrix eigvec_mat:
    arma::mat eigvec_mat = arma::mat(N,0);

    
    
    //Normalise the matrix eigec_i:
    for (int i=0; i<=N-1; i++){
        arma::vec eigvec_i = arma::vec(N).fill(0.);
        for (int j=0; j<=N-1; j++){
            eigvec_i(j)=eigvec(i,j);
      }
        arma::vec eigvec_i_n = arma::normalise(eigvec_i,1);
        eigvec_mat.insert_cols(0, eigvec_i_n);
    }
    

    
    //Analytical solution:
    //Define vectors eigval_an and eigvec_an:
    arma::vec eigval_an=arma::vec(N).fill(0.);
    arma::mat eigvec_an_mat = arma::mat(N,0);
    //Add the values of the vectors:
    for (int i=0; i<=N-1; i++){
        eigval_an(i)=(d+2*a*std::cos((i+1)*arma::datum::pi/n));
        
        arma::vec eigvec_an_i = arma::vec(N).fill(0.);
        
        for (int j=0; j<=N-1; j++){
            eigvec_an_i(j)=std::sin((i+1)*(j+1)*arma::datum::pi/n);
            
      }
        arma::vec eigvec_an_i_n = arma::normalise(eigvec_an_i,1);
        eigvec_an_mat.insert_cols(0, eigvec_an_i_n);
    }
    //Normalise the vector eigval_an:
    arma::vec eigval_an_n = arma::normalise(eigval_an,1);
    
    
    
    
    // Create a matrix with all the results:
    arma::mat results = arma::mat(N, 2+2*N, arma::fill::ones);
    for (int i=0; i<=N-1; i++){
        results(i,0) = eigval_n(i);
        results(i,1) = eigval_an_n(i);
        for (int j=0; j<=N-1; j++){
            results(i,2+2*j) = eigvec_mat(i,j);
            results(i,3+2*j) = eigvec_an_mat(i,j);
        }
    }
    
    std::cout << results;
        
    return 0;
}
