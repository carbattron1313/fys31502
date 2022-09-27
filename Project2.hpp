//
//  Project2.hpp
//  prueba1
//
//  Created by Elena Muñoz Rivas on 17/9/22.
//

#ifndef Project2_hpp
#define Project2_hpp
#include <armadillo>
#include <stdio.h>
#include <math.h>

// Create a tridiagonal matrix tridiag(a,d,e) of size n*n,
// from scalar input a, d, and e. That is, create a matrix where
// - all n-1 elements on the subdiagonal have value a
// - all n elements on the diagonal have value d
// - all n-1 elements on the superdiagonal have value e
arma::mat create_tridiagonal(int N, double a, double d, double e)
{
  // Start from identity matrix
  arma::mat A = arma::mat(N, N, arma::fill::eye);

  // Fill the first row (row index 0), e.g.
  A(0,0) = d;
  A(0,1) = e;

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
    for (int i=1; i<=N-2; i++){
        A(i,i+1) = e;
        A(i,i-1) = a;
        A(i,i) = d;
    }
    
  // Fill last row (row index n-1)
    A(N-1,N-1) = d;
    A(N-1,N-2) = a;
  return A;
}


// A function that finds the max off-diag element of a symmetric matrix A.
// - The matrix indices of the max element are returned by writing to the
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
  // Get size of the matrix A:
    int sizeA = A.n_rows;
  // Check that A is square and larger than 1x1. If the condition evaluates to false, the program is killed with
    // an assertion error.
    
    assert(A.is_square()==true);
    assert(sizeA>1);

  // Initialize a double variable 'maxval' to A(k,l). We'll use this variable to keep track of the largest off-diag element.
    double maxval;
  // Loop through all elements in the upper triangle of A (not including the diagonal)
    for (int i=0; i<=sizeA-2; i++){
        for (int j=1+i; j<=sizeA-1; j++){
            if(std::abs(A(i,j))> std::abs(maxval)){
                k=i;
                l=j;
                maxval = A(i,j);
            }
        }
    }
  // When encountering a matrix element with larger absolute value than the current value of maxval,
  // update k, l and max accordingly.

    return maxval;
}



void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int N){
    //Define variables s and c (sin and cos):
    double s,c;
    //Compute tan, cos, sin:
    if (A(k,l) != 0.){
        //Define varibles tau and t (tan):
        double tau, t;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0){
            t = 1./(tau + std::sqrt(1+pow(tau,2)));
        }
        else{
            t = -1./(-tau + std::sqrt(1+pow(tau,2)));
        }
        c = 1./std::sqrt(1+pow(t,2));
        s = c*t;
    } else {
        c = 1.;
        s = 0.;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    //Transform current A matrix, changing the matrix elements with indices l and k:
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0; // hard-coding of the zeros
    A(l,k) = 0.0;
    //Change the remaining elements:
    for ( int i = 0; i < N; i++ ) {
        if ( i != k && i != l ) {
        a_ik = A(i,k);
        a_il = A(i,l);
        A(i,k) = c*a_ik - s*a_il;
        A(k,i) = A(i,k);
        A(i,l) = c*a_il + s*a_ik;
        A(l,i) = A(i,l);
        }
        //Update the overall rotation matrix:
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}



//Jacobi method eigensolver:
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        int maxiter, int& iterations, bool& converged){
    //Variable iterations, counts the number of iterations:
    iterations = 0;
    //Bool that registers if the program has reached the maximun number of iterations:
    converged = true;
    //Define variables, k, l and N:
    int k;
    int l;
    int N = A.n_rows;
    
    //Define matrice R:
    arma::mat R = arma::mat(N,N, arma::fill::eye);
    
    //It Runs jacobo_rotate until max off-diagonal element < eps:
    while (std::abs(max_offdiag_symmetric(A, k, l))>eps and iterations<= maxiter){
        iterations += 1;
        max_offdiag_symmetric(A, k, l);
        jacobi_rotate(A, R, k, l, N);
    }
    for (int i=0; i<N; i++){
        //Writes the eigenvalues as entries in the vector "eigenvalues" and
        //the eigenvectors as columns in the matrix "eigenvectors"
        eigenvalues(i) = A(i,i);
        for (int j=0; j<N; j++){
            eigenvectors(i,j)=R(j,i);
        }
    }
    if (iterations<= maxiter){
        converged = false;
    }
    return;
}


#endif /* Project2_hpp */
