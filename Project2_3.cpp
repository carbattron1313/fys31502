//
//  Project2_3.cpp
//  
//
//  Created by Elena Muñoz Rivas on 21/9/22.
//

#include "Project2.hpp"

int main(){
    //Identify the largest off-diagonal element of a matrix A:
    arma::mat A = arma::mat(4, 4, arma::fill::eye);
    A(3,0) = 0.5;
    A(0,3) = 0.5;
    A(1,2) = -0.7;
    A(2,1) = -0.7;
    int k= 0;
    int l= 0;
    
    std:: cout << max_offdiag_symmetric(A, k, l) <<std::endl;
    return 0;
}

