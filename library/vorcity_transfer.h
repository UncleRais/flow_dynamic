#pragma once

#include <algorithm>

#include "utils.h"
#include "slae_solver.h"


namespace vorcity_transfer {


    void set_velocity_on_boundary(const problem_params &ps, Vector &velocity_u, Vector &velocity_v) {
        // iterations over all boundary grid nodes
        const uint Nx = ps._Nx;
        const uint Ny = ps._Ny;

        // - Bottom boundary -
        for (uint i = 0; i <= Nx; ++i) {
            const auto k = ps.ij_to_k(i, 0);
            velocity_u[k] = ps._boundary_velocities[0];
        }

        // - Right boundary -
        for (uint j = 0; j <= Ny; ++j) { 
            const auto k = ps.ij_to_k(Nx, j);
            velocity_v[k] = ps._boundary_velocities[1];
        }

        // - Top boundary -
        for (uint i = 0; i <= Nx; ++i) {
            const auto k = ps.ij_to_k(i, Ny);
            velocity_u[k] = ps._boundary_velocities[2];
        }

        // - Left boundary -
        for (uint j = 0; j <= Ny; ++j) {
            const auto k = ps.ij_to_k(0, j);
            velocity_v[k] = ps._boundary_velocities[3];
        }    
    }


    void calc_velocity_on_vertices(const problem_params& ps, const Vector& sol_prev, Vector& velocity_u, Vector& velocity_v) {
        // --- Internal vertices ---
        const uint Nx = ps._Nx;
        const uint Ny = ps._Ny;

        const T half_of_div_hx = 0.5 / ps._hx;
        const T half_of_div_hy = 0.5 / ps._hy;

        for (uint i = 1; i <= Nx - 1; ++i) {
            for (uint j = 1; j <= Ny - 1; ++j) {
                const uint k = ps.ij_to_k(i, j);

                const T psi_iplus_j  = sol_prev[ps.ij_to_row(i + 1, j,     Equation::SECOND)];
                const T psi_iminus_j = sol_prev[ps.ij_to_row(i - 1, j,     Equation::SECOND)];
                const T psi_i_jplus  = sol_prev[ps.ij_to_row(i,     j + 1, Equation::SECOND)];
                const T psi_i_jminus = sol_prev[ps.ij_to_row(i,     j - 1, Equation::SECOND)];

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
    void build_system_transfer(
        const problem_params &ps,
        std::vector<Triplet> &coef_triplets,
        const Vector &sol_prev,
        const Vector &velocity_u, const Vector &velocity_v,
        Vector &rhs
    ) {
        const auto equation = Equation::FIRST; // => filling the top half of the matrix

        const uint Nx      = ps._Nx;
        const uint Ny      = ps._Ny;
        const uint count_x = ps._count_x;
        const uint count_y = ps._count_y;

        const T div_tau = 1.0 / ps._tau;
        const T div_hx  = 1.0 / ps._hx;
        const T div_hy  = 1.0 / ps._hy;

        const T div_hx2 = div_hx * div_hx;
        const T div_hy2 = div_hy * div_hy;

        const T nu_div_hx2 = ps._nu * div_hx2;
        const T nu_div_hy2 = ps._nu * div_hy2;

        // --- Internal vertices ---
        // Iterate regular 'cross-shaped' template

        for (uint i = 1; i <= Nx - 1; ++i) {
            for (uint j = 1; j <= Ny - 1; ++j) {
                // Fill RHS
                const auto row = ps.ij_to_row(i, j, equation);

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
                const T a_S =  0.5 * div_hy * velocity_v[ps.ij_to_k(i,     j - 1)] - nu_div_hy2;
                const T a_N = -0.5 * div_hy * velocity_v[ps.ij_to_k(i,     j + 1)] - nu_div_hy2;

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
        for (uint i = 0; i <= Nx; ++i) {
            // Fill RHS
            const auto row = ps.ij_to_row(i, 0, equation);

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
        for (uint j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            const auto row = ps.ij_to_row(Nx, j, equation);

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
        for (uint i = 0; i <= Nx; ++i) {
            // Fill RHS
            const auto row = ps.ij_to_row(i, Ny, equation);

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
        for (uint j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            const auto row = ps.ij_to_row(0, j, equation);

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
    void build_system_laplace(
        const problem_params &ps,
        std::vector<Triplet> &coef_triplets,
        const Vector &sol_prev,
        const Vector &velocity_u, const Vector &velocity_v,
        Vector &rhs
    ) {
        const auto equation = Equation::SECOND; // => filling the bottom half of the matrix

        const uint Nx      = ps._Nx;
        const uint Ny      = ps._Ny;
        const uint count_x = ps._count_x;
        const uint count_y = ps._count_y;

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

        for (uint i = 1; i <= Nx - 1; ++i) {
            for (uint j = 1; j <= Ny - 1; ++j) {
                // Fill RHS
                const auto row = ps.ij_to_row(i, j, equation);

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
        for (uint i = 0; i <= Nx; ++i) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_i_0 = ps.ij_to_rc(i, 0, i, 0, equation, Variable::PSI);

            coef_triplets.push_back(Triplet(idx_psi_i_0.row, idx_psi_i_0.col, 1.0));
        }

        // - Right boundary (no corners) -
        for (uint j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_Nx_j = ps.ij_to_rc(Nx, j, Nx, j, equation, Variable::PSI);

            coef_triplets.push_back(Triplet(idx_psi_Nx_j.row, idx_psi_Nx_j.col, 1.0));
        }

        // - Top boundary -
        for (uint i = 0; i <= Nx; ++i) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_i_Ny = ps.ij_to_rc(i, Ny, i, Ny, equation, Variable::PSI);

            coef_triplets.push_back(Triplet(idx_psi_i_Ny.row, idx_psi_i_Ny.col, 1.0));
        }

        // - Left boundary (no corners) -
        for (uint j = 1; j <= Ny - 1; ++j) {
            // Fill RHS
            // < always stays zero, no need to fill every time >

            // Fill Matrix
            const auto idx_psi_0_j = ps.ij_to_rc(0, j, 0, j, equation, Variable::PSI);
            coef_triplets.push_back(Triplet(idx_psi_0_j.row, idx_psi_0_j.col, 1.0));
        }
    }


    void solve(
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
        slae_solvers::sparse_LU(coef_triplets, sol_curr, rhs, verbose);
        sol_curr.swap(sol_prev);
    }


    inline T integrate_value_over_G(const problem_params &ps, const Vector &values, std::function<T(T)> func) {
        const uint Nx = ps._Nx;
        const uint Ny = ps._Ny;

        const uint size = ps._size;

        // Compure epsilon = 0.5 * integral[u^2 + v^2, G]
        T sum = 0.0;

        for (uint i = 0; i <= Nx - 1; ++i) {
            for (uint j = 0; j <= Ny - 1; ++j) {
                const auto idx_bl = ps.ij_to_row(i, j, Equation::FIRST);
                const auto idx_br = ps.ij_to_row(i + 1, j, Equation::FIRST);
                const auto idx_tr = ps.ij_to_row(i + 1, j + 1, Equation::FIRST);
                const auto idx_tl = ps.ij_to_row(i, j + 1, Equation::FIRST);

                const T val_bl = values[idx_bl];
                const T val_br = values[idx_br];
                const T val_tr = values[idx_tr];
                const T val_tl = values[idx_tl];

                const T area = ps._hx * ps._hy;

                sum += 0.5 * (func(val_bl) + func(val_br) + func(val_tr) + func(val_tl)) * area * 0.25;
            }
        }

        return sum;
    }


    T get_integral_omega(const problem_params &ps, const Vector &sol) {
        const uint size = ps._size;

        // Extract omega from solution
        Vector omega(size / 2);
        std::copy(sol.begin(), sol.begin() + size / 2, omega.begin());

        const auto identity = [](T x) -> T { return x; };

        // Compute sum = integral[omega, G]
        const T sum = integrate_value_over_G(ps, omega, identity);

        // Normalize over max abs(psi)
        /*std::sort(omega.begin(), omega.end());
        const T max = *(omega.end() - 1);

        return (max > 0.0) ? sum / max : 0.0;*/

        return sum;
    }


    T get_integral_omega2(const problem_params &ps, const Vector &sol) {
        const uint size = ps._size;

        // Extract omega from solution
        Vector omega(size / 2);
        std::copy(sol.begin(), sol.begin() + size / 2, omega.begin());

        const auto sqr = [](T x) -> T { return x * x; };

        // Compute sum = integral[omega, G]
        const T sum = integrate_value_over_G(ps, omega, sqr);

        return sum;
    }


    T get_integral_dw_dt_psi(const problem_params &ps, const Vector &sol_prev, const Vector &sol_next) {
        const uint size = ps._size;

        // Extract psi, omega from solution
        Vector psi_prev(size / 2);
        std::copy(sol_prev.begin() + size / 2, sol_prev.end(), psi_prev.begin());

        Vector psi_next(size / 2);
        std::copy(sol_prev.begin() + size / 2, sol_prev.end(), psi_next.begin());

        Vector omega_prev(size / 2);
        std::copy(sol_prev.begin(), sol_prev.begin() + size / 2, omega_prev.begin());

        Vector omega_next(size / 2);
        std::copy(sol_next.begin(), sol_next.begin() + size / 2, omega_next.begin());

        Vector dw_dt_psi(size / 2);
        for (uint i = 0; i < size / 2; ++i)
            dw_dt_psi[i] += (omega_next[i] - omega_prev[i]) / ps._tau * 0.5 * (psi_prev[i] + psi_next[i]);

        const auto identity = [](T x) -> T { return x; };

        // Compute sum = integral[omega, G]
        const T I = integrate_value_over_G(ps, dw_dt_psi, identity);

        return I;
    }


    T get_integral_epsilon(const problem_params &ps, const Vector &velocity_u, const Vector &velocity_v) {

        // Compure epsilon = 0.5 * integral[u^2 + v^2, G]
        const auto sqr = [](T x) -> T { return x * x; };

        const T epsilon = 0.5 * (integrate_value_over_G(ps, velocity_u, sqr) + integrate_value_over_G(ps, velocity_v, sqr));

        return epsilon;
    }
}