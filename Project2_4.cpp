//
//  Project2_4.cpp
//  
//
//  Created by Elena Muñoz Rivas on 21/9/22.
//

#include "Project2.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>

int main(){
    
    //Set a filename: Problem2_6
    std::string filename = "Problem2_6.txt";
    
    std::ofstream ofile;
    ofile.open(filename);
    
    //To set up the tridiagonal A for N=6, we create the constats:
    int N = 6;
    int n = N+1;
    double h = 1./n;

    //Create the values of the vectors of the diagonal and triagonal for matrix A:
    double a = -1./(h*h);
    double d = 2./(h*h);
    arma::mat A = create_tridiagonal(N,a,d,a);
    
    //Create all the variables:
    double eps = 1e-8;
    arma::vec eigenvalues = arma::vec(N).fill(0.);
    arma::mat eigenvectors = arma::mat(N,N).fill(0.);
    int maxiter = pow(N, 3);
    int iterations;
    bool converged;
    
    
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    
    std::cout << A << std::endl
            << eigenvalues << std::endl
            << eigenvectors << std::endl
            << iterations << std::endl
            << converged << std::endl;
    
    ofile << std::setw(10) << std::setprecision(3) << eigenvalues
        << std::setw(10) << std::setprecision(3) << eigenvectors << std::endl;
    
    ofile.close();
    
    
}
