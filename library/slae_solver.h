#pragma once
#include "utils.h"

namespace slae_solver {
    // coefficients - list of non-zeros coefficients
    // rhs - the right hand side-vector resulting from the constraints
    void sparse_simplicial_cholesky(const std::vector<triplet>& coefficients, vector& solution, const vector& rhs, bool verbose = false) {
        assert(rhs.rows() > 0 && rhs.cols() > 0 && coefficients.size() > 0);
        const size_type n = rhs.rows();
        matrix A(n, n);
        A.setFromTriplets(coefficients.begin(), coefficients.end());

        if (verbose) print_system(A, rhs);

        // Solving:
        Eigen::SimplicialCholesky<matrix> Cholesky_solver(A);  // performs a Cholesky factorization of A
        solution = Cholesky_solver.solve(rhs);                 // use the factorization to solve for the given right hand side
    };

    void sparse_LU(const std::vector<triplet>& coefficients, vector& solution, const vector& rhs, bool verbose = false) {
        assert(rhs.rows() > 0 && rhs.cols() > 0 && coefficients.size() > 0);
        const size_type n = rhs.rows();
        matrix A(n, n);
        A.setFromTriplets(coefficients.begin(), coefficients.end());

        if (verbose) print_system(A, rhs);

        Eigen::SparseLU<matrix, Eigen::COLAMDOrdering<int>> LU_solver;
        LU_solver.analyzePattern(A);     // Compute the ordering permutation vector from the structural pattern of A
        LU_solver.factorize(A);          // Compute the numerical factorization 
        solution = LU_solver.solve(rhs); // Use the factors to solve the linear system             
    };

    void test_example(bool verbose = false) {
        constexpr size_type n = 5;
        std::vector<triplet> coefficients;    
        vector b(n), x(n);                          
        build_test_problem(coefficients, b, n);
        sparse_simplicial_cholesky(coefficients, x, b, verbose);
        if (verbose) std::cout << "solution : " << x << std::endl;
    };
};