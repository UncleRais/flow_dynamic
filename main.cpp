
#include <iostream>
#include "library/Eigen/Sparse"
#include "library/Eigen/Dense"
#include <omp.h>
#include <vector>
#include <string>
#include <iomanip>

typedef Eigen::SparseMatrix<double> Matrix;
typedef Eigen::Triplet<double> T;

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n){
    for (std::size_t i = 0; i < n; ++i) {
        b[i] = 1.0;
        for (std::size_t j = 0; j < n; ++j) 
            if(i == j) {
                const double val = 1.0;
                coefficients.push_back(T(i, j, val));
            };
    };
};

int main(int argc, char** argv) {

    int n = 5;  // size of the image

    // Assembly:
    std::vector<T> coefficients;            // list of non-zeros coefficients
    Eigen::VectorXd b(n);                   // the right hand side-vector resulting from the constraints
    buildProblem(coefficients, b, n);

    Matrix A(n, n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    
    for (std::size_t i = 0; i < n; ++i) {
        std::cout << "{ ";
        for (std::size_t j = 0; j < n; ++j) {
            std::cout << std::setw(9) << std::setprecision(8) << A.coeffRef(i, j) << " ";
        };
        std::cout << " }" << std::endl;
    };
    

    // Solving:
    Eigen::SimplicialCholesky<Matrix> chol(A);  // performs a Cholesky factorization of A
    Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side

    std::cout << "{ ";
    for (auto& el: x) {
        std::cout << el << " ";
    };
    std::cout << " }" << std::endl;

    return 0;
};