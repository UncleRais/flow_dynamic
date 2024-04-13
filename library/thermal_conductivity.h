#pragma once

#include "utils.h"
#include "slae_solver.h"


namespace thermal_conductivity {


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

        const T hx2 = ps._hx * ps._hx;
        const T hy2 = ps._hy * ps._hy;

        const T kappa_hx2_tau = ps._kappa * hx2 * ps._tau;
        const T kappa_hy2_tau = ps._kappa * hy2 * ps._tau;

        const T hx2_hy2 = hx2 * hy2;
        const T hx_hy2_tau = ps._hx * hy2 * ps._tau;
        const T hx2_hy_tau = hx2 * ps._hy * ps._tau;

        // --- Internal vertices ---
        // Iterate regular 'cross-shaped' template

        for (size_type i = 1; i <= Nx - 1; ++i) {
            for (size_type j = 1; j <= Ny - 1; ++j) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);

                rhs[row] = hx2_hy2 * sol_prev[row];

                // Fill Matrix
                const auto idx_w_P = ps.ij_to_rc(i, j, i,     j,     equation, Variable::TEMPERATURE);
                const auto idx_w_W = ps.ij_to_rc(i, j, i - 1, j,     equation, Variable::TEMPERATURE);
                const auto idx_w_E = ps.ij_to_rc(i, j, i + 1, j,     equation, Variable::TEMPERATURE);
                const auto idx_w_S = ps.ij_to_rc(i, j, i,     j - 1, equation, Variable::TEMPERATURE);
                const auto idx_w_N = ps.ij_to_rc(i, j, i,     j + 1, equation, Variable::TEMPERATURE);

                const T U_plus  = (velocity_u[ps.ij_to_k(i + 1, j)] + velocity_u[ps.ij_to_k(i    , j)]) * 0.5;
                const T U_minus = (velocity_u[ps.ij_to_k(i    , j)] + velocity_u[ps.ij_to_k(i - 1, j)]) * 0.5;
                const T V_plus  = (velocity_v[ps.ij_to_k(i, j + 1)] + velocity_v[ps.ij_to_k(i, j    )]) * 0.5;
                const T V_minus = (velocity_v[ps.ij_to_k(i, j    )] + velocity_v[ps.ij_to_k(i, j - 1)]) * 0.5;

                const T a_P = hx2_hy2 + 0.5 * hx_hy2_tau * (U_plus - U_minus) + 0.5 * hx2_hy_tau * (V_plus - V_minus)
                              + 2.0 * (kappa_hx2_tau + kappa_hy2_tau);
                const T a_W = -0.5 * hx_hy2_tau * U_minus - kappa_hy2_tau;
                const T a_E =  0.5 * hx_hy2_tau * U_plus  - kappa_hy2_tau;
                const T a_S = -0.5 * hx2_hy_tau * V_minus - kappa_hx2_tau;
                const T a_N =  0.5 * hx2_hy_tau * V_plus  - kappa_hx2_tau;

                coef_triplets.push_back(Triplet(idx_w_P.row, idx_w_P.col, a_P));
                coef_triplets.push_back(Triplet(idx_w_W.row, idx_w_W.col, a_W));
                coef_triplets.push_back(Triplet(idx_w_E.row, idx_w_E.col, a_E));
                coef_triplets.push_back(Triplet(idx_w_S.row, idx_w_S.col, a_S));
                coef_triplets.push_back(Triplet(idx_w_N.row, idx_w_N.col, a_N));
            }
        }

        // --- Boundary vertices ---

        const auto& BC_bottom =  ps._boundary_temperature_conditions[0];
        const auto& BC_right  =  ps._boundary_temperature_conditions[1];
        const auto& BC_top    =  ps._boundary_temperature_conditions[2];
        const auto& BC_left   =  ps._boundary_temperature_conditions[3];
        const T one_div_c_rho = ps._kappa / ps._lambda;

        // - Bottom boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            if (BC_bottom.first == ThermalBoundaryType::TEMPERATURE) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, 0, equation);
                rhs[row] = BC_bottom.second;
                // Fill Matrix
                const auto idx_T_i_0 = ps.ij_to_rc(i, 0, i, 0, equation, Variable::TEMPERATURE);
                coef_triplets.push_back(Triplet(idx_T_i_0.row, idx_T_i_0.col, 1.));
            }
            if (BC_bottom.first == ThermalBoundaryType::FLUX){
                // Fill RHS
                // const auto row = ps.ij_to_rhs_row(i, 0, equation);
                // rhs[row] = BC_bottom.second * hx2_hy2 * ps._tau;
            }
        }

        // - Right boundary (no corners) -
        for (size_type j = 1; j <= Ny - 1; ++j) {
            const size_type i = Nx;
            if (BC_right.first == ThermalBoundaryType::TEMPERATURE) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);
                rhs[row] = BC_right.second;
                // Fill Matrix
                const auto idx_T_Nx_j = ps.ij_to_rc(i, j, i, j, equation, Variable::TEMPERATURE);
                coef_triplets.push_back(Triplet(idx_T_Nx_j.row, idx_T_Nx_j.col, 1.));
            }
            if (BC_right.first == ThermalBoundaryType::FLUX){
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);

                rhs[row] = BC_right.second * one_div_c_rho * hx_hy2_tau;

                // Fill Matrix
                const auto idx_w_P = ps.ij_to_rc(i, j, i,     j,     equation, Variable::TEMPERATURE);
                const auto idx_w_W = ps.ij_to_rc(i, j, i - 1, j,     equation, Variable::TEMPERATURE);
                const auto idx_w_S = ps.ij_to_rc(i, j, i,     j - 1, equation, Variable::TEMPERATURE);
                const auto idx_w_N = ps.ij_to_rc(i, j, i,     j + 1, equation, Variable::TEMPERATURE);

                const T U_plus  = (velocity_u[ps.ij_to_k(i    , j)]) * 0.5;
                const T U_minus = (velocity_u[ps.ij_to_k(i    , j)] + velocity_u[ps.ij_to_k(i - 1, j)]) * 0.5;
                const T V_plus  = (velocity_v[ps.ij_to_k(i, j + 1)] + velocity_v[ps.ij_to_k(i, j    )]) * 0.5;
                const T V_minus = (velocity_v[ps.ij_to_k(i, j    )] + velocity_v[ps.ij_to_k(i, j - 1)]) * 0.5;

                const T a_P = hx2_hy2 + 0.5 * hx_hy2_tau * (U_plus - U_minus) + 0.5 * hx2_hy_tau * (V_plus - V_minus)
                              + 2.0 * (kappa_hx2_tau) + kappa_hy2_tau;
                const T a_W = -0.5 * hx_hy2_tau * U_minus - kappa_hy2_tau;
                const T a_S = -0.5 * hx2_hy_tau * V_minus - kappa_hx2_tau;
                const T a_N =  0.5 * hx2_hy_tau * V_plus  - kappa_hx2_tau;

                coef_triplets.push_back(Triplet(idx_w_P.row, idx_w_P.col, a_P));
                coef_triplets.push_back(Triplet(idx_w_W.row, idx_w_W.col, a_W));
                coef_triplets.push_back(Triplet(idx_w_S.row, idx_w_S.col, a_S));
                coef_triplets.push_back(Triplet(idx_w_N.row, idx_w_N.col, a_N));
            }
        }

        // - Top boundary -
        for (size_type i = 0; i <= Nx; ++i) {
            if (BC_top.first == ThermalBoundaryType::TEMPERATURE) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, Ny, equation);
                rhs[row] = BC_top.second;
                // Fill Matrix
                const auto idx_T_i_Ny = ps.ij_to_rc(i, Ny, i, Ny, equation, Variable::TEMPERATURE);
                coef_triplets.push_back(Triplet(idx_T_i_Ny.row, idx_T_i_Ny.col, 1.));
            }
            if (BC_top.first == ThermalBoundaryType::FLUX){
                // Fill RHS
                // const auto row = ps.ij_to_rhs_row(i, 0, equation);
                // rhs[row] = BC_top.second * hx2_hy2 * ps._tau;
            }
        }

        // - Left boundary (no corners) -
        for (size_type j = 1; j <= Ny - 1; ++j) {
            const size_type i = 0;
            if (BC_left.first == ThermalBoundaryType::TEMPERATURE) {
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);
                rhs[row] = BC_left.second;
                // Fill Matrix
                const auto idx_T_0_j = ps.ij_to_rc(i, j, i, j, equation, Variable::TEMPERATURE);
                coef_triplets.push_back(Triplet(idx_T_0_j.row, idx_T_0_j.col, 1.));
            }
            if (BC_left.first == ThermalBoundaryType::FLUX){
                // Fill RHS
                const auto row = ps.ij_to_rhs_row(i, j, equation);

                rhs[row] = BC_left.second * one_div_c_rho * hx_hy2_tau;

                // Fill Matrix
                const auto idx_w_P = ps.ij_to_rc(i, j, i,     j,     equation, Variable::TEMPERATURE);
                const auto idx_w_E = ps.ij_to_rc(i, j, i + 1, j,     equation, Variable::TEMPERATURE);
                const auto idx_w_S = ps.ij_to_rc(i, j, i,     j - 1, equation, Variable::TEMPERATURE);
                const auto idx_w_N = ps.ij_to_rc(i, j, i,     j + 1, equation, Variable::TEMPERATURE);

                const T U_plus  = (velocity_u[ps.ij_to_k(i + 1, j)]     + velocity_u[ps.ij_to_k(i, j)]) * 0.5;
                const T U_minus = (velocity_u[ps.ij_to_k(i,     j)]) * 0.5;
                const T V_plus  = (velocity_v[ps.ij_to_k(i,     j + 1)] + velocity_v[ps.ij_to_k(i, j    )]) * 0.5;
                const T V_minus = (velocity_v[ps.ij_to_k(i,     j    )] + velocity_v[ps.ij_to_k(i, j - 1)]) * 0.5;

                const T a_P = hx2_hy2 + 0.5 * hx_hy2_tau * (U_plus - U_minus) + 0.5 * hx2_hy_tau * (V_plus - V_minus)
                              + 2.0 * (kappa_hx2_tau) + kappa_hy2_tau;
                const T a_E =  0.5 * hx_hy2_tau * U_plus  - kappa_hy2_tau;
                const T a_S = -0.5 * hx2_hy_tau * V_minus - kappa_hx2_tau;
                const T a_N =  0.5 * hx2_hy_tau * V_plus  - kappa_hx2_tau;

                coef_triplets.push_back(Triplet(idx_w_P.row, idx_w_P.col, a_P));
                coef_triplets.push_back(Triplet(idx_w_E.row, idx_w_E.col, a_E));
                coef_triplets.push_back(Triplet(idx_w_S.row, idx_w_S.col, a_S));
                coef_triplets.push_back(Triplet(idx_w_N.row, idx_w_N.col, a_N));
            }    
        }
    }


    inline void solve (
        const problem_params &ps,
        Vector &sol_prev,   Vector &sol_curr,
        const Vector &velocity_u, const Vector &velocity_v,
        Vector &rhs
    ) {

        std::vector<Triplet> coef_triplets; // non-zero coefs of matrix 'A' as { i, j, value } triplets
        coef_triplets.reserve(1 * ps._size); // not an accurate estimate

        build_system_transfer(ps, coef_triplets, sol_prev, velocity_u, velocity_v, rhs);

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