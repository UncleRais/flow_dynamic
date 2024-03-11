#pragma once



#include "utils.h"


namespace slae {

    inline void validate_system(const Matrix &A, const Vector &solution, const Vector &rhs) {
        const auto size = rhs.rows();

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

    inline void sparse_solve_decomposition(
        const std::vector<Triplet> &coef_triplets,
        Vector &solution,
        const Vector &rhs,
        DecompositionMethod method,
        bool verbose = false
    ) {
        const auto size = rhs.rows();
        
        // Build system matrix
        Matrix A(size, size);
        A.setFromTriplets(coef_triplets.begin(), coef_triplets.end());
        A.makeCompressed();

        // Check system validity
        if (verbose) std::ofstream("log.txt") << system_as_string(A, rhs);
        slae::validate_system(A, solution, rhs);

        // Solve
        if (method == DecompositionMethod::SPARSE_LU) {
            Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>> LU_solver;

            LU_solver.analyzePattern(A);     // compute the ordering permutation vector from the structural pattern of A
            LU_solver.factorize(A);          // compute the numerical factorization 
            
            solution = LU_solver.solve(rhs); // use factorization to solve for rhs      
        }
        if (method == DecompositionMethod::SPARSE_QR) {
            Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> QR_solver;

            QR_solver.analyzePattern(A);
            QR_solver.factorize(A);
           
            solution = QR_solver.solve(rhs);
        }
        if (method == DecompositionMethod::SIMPLICIAL_CHOLESKY) {
            Eigen::SimplicialCholesky<Matrix> Cholesky_solver(A);  // performs a Cholesky factorization of A

            Cholesky_solver.analyzePattern(A);
            Cholesky_solver.factorize(A);
            
            solution = Cholesky_solver.solve(rhs);                 // uses the factorization to solve for rhs
        }
    }

    inline void sparse_solve_iterative(
        const std::vector<Triplet> &coef_triplets,
        Vector &solution,
        const Vector &rhs,
        IterativeMethodData method_data,
        Vector* initial_guess = nullptr, 
        bool verbose = false
    ) {
        const auto size = rhs.rows();

        // Build system matrix
        Matrix A(size, size);
        A.setFromTriplets(coef_triplets.begin(), coef_triplets.end());
        A.makeCompressed();

        // Check system validity
        if (verbose) std::ofstream("log.txt") << system_as_string(A, rhs);
        slae::validate_system(A, solution, rhs);

        // Solve
        if (method_data.method == IterativeMethod::LS_CONJUGATE_GRADIENT) {
            Eigen::LeastSquaresConjugateGradient<Matrix> ConjGrad_solver;

            if (method_data.tolerance.has_value()) ConjGrad_solver.setTolerance(method_data.tolerance.value());
            if (method_data.max_iterations.has_value()) ConjGrad_solver.setMaxIterations(method_data.max_iterations.value());

            ConjGrad_solver.compute(A);

            if (initial_guess)
                solution = ConjGrad_solver.solveWithGuess(rhs, *initial_guess);
            else
                solution = ConjGrad_solver.solve(rhs);
        }
        if (method_data.method == IterativeMethod::BI_CGSTAB) {
            Eigen::BiCGSTAB<Matrix> BiCGSTAB_solver;

            if (method_data.tolerance.has_value()) BiCGSTAB_solver.setTolerance(method_data.tolerance.value());
            if (method_data.max_iterations.has_value()) BiCGSTAB_solver.setMaxIterations(method_data.max_iterations.value());
            
            BiCGSTAB_solver.compute(A);
            
            if (initial_guess)
                solution = BiCGSTAB_solver.solveWithGuess(rhs, *initial_guess);
            else
                solution = BiCGSTAB_solver.solve(rhs);
        }
    }
}