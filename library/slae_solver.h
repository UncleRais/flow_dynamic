#pragma once

#include "utils.h"


namespace slae_solvers {


    inline void validate_system(const Matrix &A, const Vector &solution, const Vector &rhs) {
        const uint size = rhs.rows();

        const auto get_sparse_matrix_rank = [](const Matrix &matrix) -> Eigen::Index {
            Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> QR_solver;
            QR_solver.analyzePattern(matrix);
            QR_solver.factorize(matrix);
            const auto rank = QR_solver.rank();
            return rank;
        }; // doesn't get evaluated in 'Release'
        
        assert(size == A.rows());
        assert(size == A.cols());
        assert(size == rhs.rows());
        assert(size == solution.rows());
        assert(size == get_sparse_matrix_rank(A));
    }

    void sparse_simplicial_cholesky(const std::vector<Triplet> &coef_triplets, Vector &solution, const Vector &rhs, bool verbose = false) {
        const uint size = rhs.rows();

        // Build system matrix
        Matrix A(size, size);
        A.setFromTriplets(coef_triplets.begin(), coef_triplets.end());
        A.makeCompressed();

        // Check system validity
        if (verbose) std::cout << system_as_string(A, rhs);
        slae_solvers::validate_system(A, solution, rhs);

        // Solve
        Eigen::SimplicialCholesky<Matrix> Cholesky_solver(A);  // performs a Cholesky factorization of A
        solution = Cholesky_solver.solve(rhs);                 // uses the factorization to solve for rhs
    }


    void sparse_LU(const std::vector<Triplet> &coef_triplets, Vector &solution, const Vector &rhs, bool verbose = false) {
        const uint size = rhs.rows();
        
        // Build system matrix
        Matrix A(size, size);
        A.setFromTriplets(coef_triplets.begin(), coef_triplets.end());
        A.makeCompressed();

        // Check system validity
        if (verbose) std::ofstream("log.txt") << system_as_string(A, rhs);
        slae_solvers::validate_system(A, solution, rhs);

        // Solve
        Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>> LU_solver;

        LU_solver.analyzePattern(A);     // computes the ordering permutation vector from the structural pattern of A
        LU_solver.factorize(A);          // computes the numerical factorization 
        solution = LU_solver.solve(rhs); // uses factorization to solve for rhs         
    }
}