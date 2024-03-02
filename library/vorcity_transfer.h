#pragma once
#include "utils.h"
#include "slae_solver.h"

namespace vorcity_transfer {
    using namespace slae_solver;

    void build_system(const problem_params& ps, std::vector<triplet>& A_nonzero, vector& rhs) {
        // iterations over all internal grid nodes -> regular scheme template
        const size_type Nx = ps._Nx, Ny = ps._Ny; 
        for (size_type i = 1; i < Nx; ++i) {
            for (size_type j = 1; i < Ny; ++i) {
                //auto idx = ps.omega_ij_to_rowcol(i, j, i, j, 0);
                //A_nonzero.push_back()

            };
        };

        // iterations over all boundary grid nodes -> boundary condition

    };

    void calc_velocity_on_edges() {

    };

    void calc_velocity_on_vertices() {

    };

    void solve(const problem_params& ps, vector& sol_prev, vector& sol_curr, vector& rhs) {
        std::vector<triplet> A_nonzero;                             
        build_system(ps, A_nonzero, rhs);
        sparse_simplicial_cholesky(A_nonzero, sol_curr, rhs, false);
        sol_curr.swap(sol_prev);
    };
};