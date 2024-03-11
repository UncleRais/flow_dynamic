#pragma once

#include <algorithm>

#include "utils.h"
#include "slae_solver.h"


namespace vorcity_transfer {


    inline void set_velocity_on_boundary(const problem_params &ps, Vector &velocity_u, Vector &velocity_v) {
        // iterations over all boundary grid nodes
        const size_type Nx = ps._Nx;
        const size_type Ny = ps._Ny;

        // - Bottom boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            const auto k = ps.ij_to_k(i, 0);
            velocity_u[k] = ps._boundary_velocities[0];
            //velocity_v[k] = T(0.0);
        }

        // - Right boundary -
        for (size_type j = 0; j <= Ny; ++j) { 
            const auto k = ps.ij_to_k(Nx, j);
            //velocity_u[k] = T(0.0);
            velocity_v[k] = ps._boundary_velocities[1];
        }

        // - Top boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            const auto k = ps.ij_to_k(i, Ny);
            velocity_u[k] = ps._boundary_velocities[2];
            //velocity_v[k] = T(0.0);
        }

        // - Left boundary -
        for (size_type j = 0; j <= Ny; ++j) {
            const auto k = ps.ij_to_k(0, j);
            //velocity_u[k] = T(0.0);
            velocity_v[k] = ps._boundary_velocities[3];
        }    
    }


    inline void calc_velocity_on_vertices(const problem_params& ps, const Vector& sol_prev, Vector& velocity_u, Vector& velocity_v) {
        // --- Internal vertices ---
        const size_type Nx = ps._Nx;
        const size_type Ny = ps._Ny;

        const T half_of_div_hx = 0.5 / ps._hx;
        const T half_of_div_hy = 0.5 / ps._hy;

        for (size_type i = 1; i <= Nx - 1; ++i) {
            for (size_type j = 1; j <= Ny - 1; ++j) {
                const size_type k = ps.ij_to_k(i, j);

                const T psi_iplus_j  = sol_prev[ps.ij_to_sol_row(i + 1, j,     Variable::PSI)];
                const T psi_iminus_j = sol_prev[ps.ij_to_sol_row(i - 1, j,     Variable::PSI)];
                const T psi_i_jplus  = sol_prev[ps.ij_to_sol_row(i,     j + 1, Variable::PSI)];
                const T psi_i_jminus = sol_prev[ps.ij_to_sol_row(i,     j - 1, Variable::PSI)];

                velocity_u[k] =  half_of_div_hy * (psi_i_jplus - psi_i_jminus);
                velocity_v[k] = -half_of_div_hx * (psi_iplus_j - psi_iminus_j);
            }
        }
    }


    // --- 1st equation ---
    //
    // Template on internal points: 'crosss-shaped'
    //
    //      S                i_j+1
    //      |                 |
    //  W---P---E    i-1_j---i_j---i+1_j
    //      |                 |
    //      N                i_j-1
    //
    // Geographic notation:
    //
    // P ~ i_j
    // W ~ i-1_j
    // E ~ i+1_j
    // S ~ i_j-1
    // N ~ i_j+1
    //
    // Template equation:
    // a_P * w_P  +  a_W * w_W  +  a_E * w_E  +  a_S * w_S  + a_N * w_N = f_i
    //
    // While iterating we use .ij_to_rc() and .it_to_rhs_row() methods to automatically
    // compute 'row' & 'col' matrix indexation corresponding to each term of the template
    //
    inline void build_system_transfer(
        const problem_params &ps,
        std::vector<Triplet> &coef_triplets,
        const Vector &sol_prev,
        const Vector &velocity_u, const Vector &velocity_v,
        Vector &rhs
    ) {
        const auto equation = Equation::FIRST; // => filling the top half of the matrix

        const size_type Nx      = ps._Nx;
        const size_type Ny      = ps._Ny;
        const size_type count_x = ps._count_x;
        const size_type count_y = ps._count_y;

        const T div_tau = 1.0 / ps._tau;
        const T div_hx  = 1.0 / ps._hx;
        const T div_hy  = 1.0 / ps._hy;

        const T div_hx2 = div_hx * div_hx;
        const T div_hy2 = div_hy * div_hy;

        const T nu_div_hx2 = ps._nu * div_hx2;
        const T nu_div_hy2 = ps._nu * div_hy2;

        // --- Internal vertices ---
        // Iterate regular 'cross-shaped' template

        for (size_type i = 1; i <= Nx - 1; ++i) {
            for (size_type j = 1; j <= Ny - 1; ++j) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);

                rhs[row] = div_tau * sol_prev[row]; /// Temperature impact will be added here

                // Fill Matrix
                const auto idx_w_P = ps.ij_to_rc(i, j, i,     j,     equation, Variable::OMEGA);
                const auto idx_w_W = ps.ij_to_rc(i, j, i - 1, j,     equation, Variable::OMEGA);
                const auto idx_w_E = ps.ij_to_rc(i, j, i + 1, j,     equation, Variable::OMEGA);
                const auto idx_w_S = ps.ij_to_rc(i, j, i,     j - 1, equation, Variable::OMEGA);
                const auto idx_w_N = ps.ij_to_rc(i, j, i,     j + 1, equation, Variable::OMEGA);

                const T a_P = div_tau + 2.0 * (nu_div_hx2 + nu_div_hy2);
                const T a_W = -0.5 * div_hx * velocity_u[ps.ij_to_k(i - 1, j    )] - nu_div_hx2;
                const T a_E =  0.5 * div_hx * velocity_u[ps.ij_to_k(i + 1, j    )] - nu_div_hx2;
                const T a_S = -0.5 * div_hy * velocity_v[ps.ij_to_k(i,     j - 1)] - nu_div_hy2;
                const T a_N =  0.5 * div_hy * velocity_v[ps.ij_to_k(i,     j + 1)] - nu_div_hy2;

                coef_triplets.push_back(Triplet(idx_w_P.row, idx_w_P.col, a_P));
                coef_triplets.push_back(Triplet(idx_w_W.row, idx_w_W.col, a_W));
                coef_triplets.push_back(Triplet(idx_w_E.row, idx_w_E.col, a_E));
                coef_triplets.push_back(Triplet(idx_w_S.row, idx_w_S.col, a_S));
                coef_triplets.push_back(Triplet(idx_w_N.row, idx_w_N.col, a_N));
            }
        }

        // --- Boundary vertices ---
        // Apply Thom conditions

        const T two_div_hx_sq = 2.0 * div_hx2;
        const T two_div_hy_sq = 2.0 * div_hy2;

        const T BVI_bottom =  2.0 * div_hx * ps._boundary_velocities[0];
        const T BVI_right  =  2.0 * div_hy * ps._boundary_velocities[1];
        const T BVI_top    = -2.0 * div_hx * ps._boundary_velocities[2];
        const T BVI_left   = -2.0 * div_hy * ps._boundary_velocities[3];
            // 'BVI' means 'Boundary Velocity Impact', refers to a Thom condition term that includes boundary velocity

        // - Bottom boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            // Fill RHS
            const auto row = ps.ij_to_rhs_row(i, 0, equation);

            rhs[row] = BVI_bottom;

            // Fill Matrix
            const auto idx_w_i_0   = ps.ij_to_rc(i, 0, i, 0, equation, Variable::OMEGA);
            const auto idx_psi_i_0 = ps.ij_to_rc(i, 0, i, 0, equation, Variable::PSI  );
            const auto idx_psi_i_1 = ps.ij_to_rc(i, 0, i, 1, equation, Variable::PSI  );

            coef_triplets.push_back(Triplet(idx_w_i_0.row,   idx_w_i_0.col,    1.0          ));
            coef_triplets.push_back(Triplet(idx_psi_i_0.row, idx_psi_i_0.col, -two_div_hy_sq));
            coef_triplets.push_back(Triplet(idx_psi_i_1.row, idx_psi_i_1.col, +two_div_hy_sq));
        }

        // - Right boundary (no corners) -
        for (size_type j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            const auto row = ps.ij_to_rhs_row(Nx, j, equation);

            rhs[row] = BVI_right;

            // Fill Matrix
            const auto idx_w_Nx_j        = ps.ij_to_rc(Nx, j, Nx,     j, equation, Variable::OMEGA);
            const auto idx_psi_Nx_j      = ps.ij_to_rc(Nx, j, Nx,     j, equation, Variable::PSI  );
            const auto idx_psi_Nxminus_j = ps.ij_to_rc(Nx, j, Nx - 1, j, equation, Variable::PSI  );

            coef_triplets.push_back(Triplet(idx_w_Nx_j.row,        idx_w_Nx_j.col,         1.0          ));
            coef_triplets.push_back(Triplet(idx_psi_Nx_j.row,      idx_psi_Nx_j.col,      -two_div_hx_sq));
            coef_triplets.push_back(Triplet(idx_psi_Nxminus_j.row, idx_psi_Nxminus_j.col, +two_div_hx_sq));
        }

        // - Top boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            // Fill RHS
            const auto row = ps.ij_to_rhs_row(i, Ny, equation);

            rhs[row] = BVI_top; 

            // Fill Matrix
            const auto idx_w_i_Ny        = ps.ij_to_rc(i, Ny, i, Ny,     equation, Variable::OMEGA);
            const auto idx_psi_i_Ny      = ps.ij_to_rc(i, Ny, i, Ny,     equation, Variable::PSI  );
            const auto idx_psi_i_Nyminus = ps.ij_to_rc(i, Ny, i, Ny - 1, equation, Variable::PSI  );
            
            coef_triplets.push_back(Triplet(idx_w_i_Ny.row,        idx_w_i_Ny.col,         1.0          ));
            coef_triplets.push_back(Triplet(idx_psi_i_Ny.row,      idx_psi_i_Ny.col,      -two_div_hy_sq));
            coef_triplets.push_back(Triplet(idx_psi_i_Nyminus.row, idx_psi_i_Nyminus.col, +two_div_hy_sq));
        }

        // - Left boundary (no corners) -
        for (size_type j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            const auto row = ps.ij_to_rhs_row(0, j, equation);

            rhs[row] = BVI_left;

            // Fill Matrix
            const auto idx_w_0_j   = ps.ij_to_rc(0, j, 0, j, equation, Variable::OMEGA);
            const auto idx_psi_0_j = ps.ij_to_rc(0, j, 0, j, equation, Variable::PSI  );
            const auto idx_psi_1_j = ps.ij_to_rc(0, j, 1, j, equation, Variable::PSI  );

            coef_triplets.push_back(Triplet(idx_w_0_j.row,   idx_w_0_j.col,    1.0          ));
            coef_triplets.push_back(Triplet(idx_psi_0_j.row, idx_psi_0_j.col, -two_div_hx_sq));
            coef_triplets.push_back(Triplet(idx_psi_1_j.row, idx_psi_1_j.col, +two_div_hx_sq));
            
        }
    }


    // --- 2nd equation ---
    //
    // Template on internal points: 'crosss-shaped'
    //
    //      S                i_j+1
    //      |                 |
    //  W---P---E    i-1_j---i_j---i+1_j
    //      |                 |
    //      N                i_j-1
    //
    // Template equation:
    // w_P  +  a_P * psi_P  +  a_W * psi_W  +  a_E * psi_E  +  a_S * psi_S  + a_N * psi_N = 0
    //
    inline void build_system_laplace(
        const problem_params &ps,
        std::vector<Triplet> &coef_triplets,
        const Vector &sol_prev,
        const Vector &velocity_u, const Vector &velocity_v,
        Vector &rhs
    ) {
        const auto equation = Equation::SECOND; // => filling the bottom half of the matrix

        const size_type Nx      = ps._Nx;
        const size_type Ny      = ps._Ny;
        const size_type count_x = ps._count_x;
        const size_type count_y = ps._count_y;

        const T div_hx = 1.0 / ps._hx;
        const T div_hy = 1.0 / ps._hy;

        const T div_hx2 = div_hx * div_hx;
        const T div_hy2 = div_hy * div_hy;
       
        const T a_P = -2.0 * (div_hx2 + div_hy2);
        const T a_W = div_hx2;
        const T a_E = div_hx2;
        const T a_S = div_hy2;
        const T a_N = div_hy2;

        // --- Internal vertices ---
        // Iterate regular 'cross-shaped' template

        for (size_type i = 1; i <= Nx - 1; ++i) {
            for (size_type j = 1; j <= Ny - 1; ++j) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);

                rhs[row] = 0.0;

                // Fill Matrix
                const auto idx_w_P   = ps.ij_to_rc(i, j, i,     j,     equation, Variable::OMEGA);
                const auto idx_psi_P = ps.ij_to_rc(i, j, i,     j,     equation, Variable::PSI  );
                const auto idx_psi_W = ps.ij_to_rc(i, j, i - 1, j,     equation, Variable::PSI  );
                const auto idx_psi_E = ps.ij_to_rc(i, j, i + 1, j,     equation, Variable::PSI  );
                const auto idx_psi_S = ps.ij_to_rc(i, j, i,     j - 1, equation, Variable::PSI  );
                const auto idx_psi_N = ps.ij_to_rc(i, j, i,     j + 1, equation, Variable::PSI  );

                coef_triplets.push_back(Triplet(idx_w_P.row,   idx_w_P.col,   1.0));
                coef_triplets.push_back(Triplet(idx_psi_P.row, idx_psi_P.col, a_P));
                coef_triplets.push_back(Triplet(idx_psi_W.row, idx_psi_W.col, a_W));
                coef_triplets.push_back(Triplet(idx_psi_E.row, idx_psi_E.col, a_E));
                coef_triplets.push_back(Triplet(idx_psi_S.row, idx_psi_S.col, a_S));
                coef_triplets.push_back(Triplet(idx_psi_N.row, idx_psi_N.col, a_N));
            }
        }

        // --- Boundary vertices ---
        // Apply zero-boundary condition: { psi  |boundary  = 0
        
        // - Bottom boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_i_0 = ps.ij_to_rc(i, 0, i, 0, equation, Variable::PSI);

            coef_triplets.push_back(Triplet(idx_psi_i_0.row, idx_psi_i_0.col, 1.0));
        }

        // - Right boundary (no corners) -
        for (size_type j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_Nx_j = ps.ij_to_rc(Nx, j, Nx, j, equation, Variable::PSI);

            coef_triplets.push_back(Triplet(idx_psi_Nx_j.row, idx_psi_Nx_j.col, 1.0));
        }

        // - Top boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_i_Ny = ps.ij_to_rc(i, Ny, i, Ny, equation, Variable::PSI);

            coef_triplets.push_back(Triplet(idx_psi_i_Ny.row, idx_psi_i_Ny.col, 1.0));
        }

        // - Left boundary (no corners) -
        for (size_type j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_0_j = ps.ij_to_rc(0, j, 0, j, equation, Variable::PSI);
            coef_triplets.push_back(Triplet(idx_psi_0_j.row, idx_psi_0_j.col, 1.0));
        }
    }


    inline void solve (
        const problem_params &ps,
        Vector &sol_prev,   Vector &sol_curr,
        Vector &velocity_u, Vector &velocity_v,
        Vector &rhs
    ) {
        calc_velocity_on_vertices(ps, sol_prev, velocity_u, velocity_v);

        std::vector<Triplet> coef_triplets; // non-zero coefs of matrix 'A' as { i, j, value } triplets
        coef_triplets.reserve(3 * ps._size); // not an accurate estimate

        build_system_transfer(ps, coef_triplets, sol_prev, velocity_u, velocity_v, rhs);
        build_system_laplace (ps, coef_triplets, sol_prev, velocity_u, velocity_v, rhs);

        const bool verbose = (ps._size < 30);
        const auto &method = ps._method;

        //const slae::Method method = slae::DecompositionMethod::SPARSE_LU;
        //const Method method = IterativeMethodData{ IterativeMethod::BI_CGSTAB, 1e-4, 200, sol_prev };

        // Reroute SLAE args to decomposition or iterative solver
        if (std::holds_alternative<DecompositionMethod>(method)) {
            const auto &decomp_method = std::get<DecompositionMethod>(method);

            slae::sparse_solve_decomposition(coef_triplets, sol_curr, rhs, decomp_method, verbose);
        }
        if (std::holds_alternative<IterativeMethodData>(method)) {
            const auto &iterative_method_data = std::get<IterativeMethodData>(method);

            slae::sparse_solve_iterative(coef_triplets, sol_curr, rhs, iterative_method_data, &sol_prev, verbose);
        }

        sol_curr.swap(sol_prev);
    }

}