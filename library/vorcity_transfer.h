#pragma once
#include "utils.h"
#include "slae_solver.h"

namespace vorcity_transfer {
    using namespace slae_solver;

    void calc_velocity_on_vertices(const problem_params& ps, const vector& sol_prev, vector& velocity_u, vector& velocity_v) {
        // iterations over all internal grid nodes
        const size_type Nx = ps._Nx, Ny = ps._Ny; 
        const double div_hx = 0.5 / ps._hx, div_hy = 0.5 / ps._hy;
        for (size_type i = 1; i < Nx; ++i)
            for (size_type j = 1; i < Ny; ++i) {
                size_type k = ps.ij_to_k(i, j, 0);
                velocity_u[k] = div_hy * (sol_prev[ps.ij_to_k(i, j + 1, 1)] - sol_prev[ps.ij_to_k(i, j - 1, 1)]);
                velocity_v[k] = div_hx * (sol_prev[ps.ij_to_k(i + 1, j, 1)] - sol_prev[ps.ij_to_k(i - 1, j, 1)]);
            };
        // we don't need to iterate over all boundary nodes as V = 0 on boundary
    };

    //      S               ij+1
    //      |                |
    //  W---p---E    i-1j---ij---i+1j
    //      |                |
    //      N               ij-1
    // template  ai    * wij   + bi    * wi+1j + ci  * wi-1j  + di    * wij+1 + ei    * wij-1 = fi
    // ai - indexation in mesh, related with central template vertice 
    // template  a{pp} * w{p}  + b{pE} * w{E}  + c{pW} * w{W} + d{pS} * w{S}  + e{pN} * w{N}  = f{p}
    // a{rc} - indexation in matrix, r - row, c - col


    void build_system_transfer(const problem_params& ps, std::vector<triplet>& A_nonzero, const vector& sol_prev,
                               const vector& velocity_u, const vector& velocity_v, vector& rhs) {
        const size_type eq = 0; //here we are working with top part of the matrix
        const size_type Nx = ps._Nx, Ny = ps._Ny; 
        const double div_tau = 1.0 / ps._tau, div_hx = 1.0 / ps._hx, div_hy = 1.0 / ps._hy;
        const double div_hx_sq = div_hx * div_hx, div_hy_sq = div_hy * div_hy;
        const double nu_hx = ps._nu * div_tau * div_hx_sq, nu_hy = ps._nu * div_tau * div_hy_sq;
        const double app = div_tau + nu_hx + nu_hy;
        // iterations over all internal grid nodes -> regular scheme template
        for (size_type i = 1; i < Nx; ++i)
            for (size_type j = 1; i < Ny; ++i) {
                //rhs
                size_type k = ps.ij_to_k(i, j, eq);
                // HERE will be the temperature impact
                rhs[k] = div_tau * sol_prev[k];

                //template
                auto id_wp = ps.ij_to_rc(i, j, i, j, eq, 0);
                A_nonzero.push_back(triplet(id_wp.first, id_wp.second, app));

                auto id_wE = ps.ij_to_rc(i, j, i + 1, j, eq, 0);
                const double apE = 0.5 * div_hx * velocity_u[ps.ij_to_k(i + 1, j, 0)] - nu_hx;
                A_nonzero.push_back(triplet(id_wE.first, id_wE.second, apE));

                auto id_wW = ps.ij_to_rc(i, j, i - 1, j, eq, 0);
                const double apW = - 0.5 * div_hx * velocity_u[ps.ij_to_k(i - 1, j, 0)] - nu_hx;
                A_nonzero.push_back(triplet(id_wW.first, id_wW.second, apW));

                auto id_wS = ps.ij_to_rc(i, j, i, j + 1, eq, 0);
                const double apS = 0.5 * div_hy * velocity_v[ps.ij_to_k(i, j + 1, 0)] - nu_hy;
                A_nonzero.push_back(triplet(id_wS.first, id_wS.second, apS));

                auto id_wN = ps.ij_to_rc(i, j, i, j - 1, eq, 0);
                const double apN = - 0.5 * div_hy * velocity_v[ps.ij_to_k(i, j - 1, 0)] - nu_hy;
                A_nonzero.push_back(triplet(id_wN.first, id_wN.second, apN));
            };

        // iterations over all boundary grid nodes -> boundary condition
        const size_type count_x = ps._count_x, count_y = ps._count_y;
        const double two_div_hx_sq = 8.0 * div_hx_sq, two_div_hy_sq = 8.0 * div_hy_sq;
        for (size_type i = 0; i < count_x; ++i) {       // BC bot
            auto id_wi0 =    ps.ij_to_rc(i, 0, i, 0, eq, 0);
            auto id_psi_i0 = ps.ij_to_rc(i, 0, i, 0, eq, 1);
            auto id_psi_i1 = ps.ij_to_rc(i, 0, i, 1, eq, 1);
            A_nonzero.push_back(triplet(id_wi0.first,    id_wi0.second,     1.0));
            A_nonzero.push_back(triplet(id_psi_i1.first, id_psi_i1.second, -two_div_hy_sq));
            A_nonzero.push_back(triplet(id_psi_i0.first, id_psi_i0.second,  two_div_hy_sq));
        };

        for (size_type j = 1; j < count_y - 1; ++j) {   // BC left (without corners)
            auto id_w0j =    ps.ij_to_rc(0, j, 0, j, eq, 0);
            auto id_psi_0j = ps.ij_to_rc(0, j, 0, j, eq, 1);
            auto id_psi_1j = ps.ij_to_rc(0, j, 1, j, eq, 1);
            A_nonzero.push_back(triplet(id_w0j.first,    id_w0j.second,     1.0));
            A_nonzero.push_back(triplet(id_psi_1j.first, id_psi_1j.second, -two_div_hx_sq));
            A_nonzero.push_back(triplet(id_psi_0j.first, id_psi_0j.second,  two_div_hx_sq));
        };

        for (size_type i = 0; i < count_x; ++i) {       // BC top
            auto id_wiNy =          ps.ij_to_rc(i, Ny, i, Ny,     eq, 0);
            auto id_psi_iNy_minus = ps.ij_to_rc(i, Ny, i, Ny - 1, eq, 1);
            auto id_psi_iNy =       ps.ij_to_rc(i, Ny, i, Ny,     eq, 1);
            A_nonzero.push_back(triplet(id_wiNy.first,          id_wiNy.second,           1.0));
            A_nonzero.push_back(triplet(id_psi_iNy.first,       id_psi_iNy.second,        two_div_hy_sq));
            A_nonzero.push_back(triplet(id_psi_iNy_minus.first, id_psi_iNy_minus.second, -two_div_hy_sq));
        };

        for (size_type j = 1; j < count_y - 1; ++j) {   // BC right (without corners)
            auto id_wNxj =          ps.ij_to_rc(Nx, j, Nx,     j, eq, 0);
            auto id_psi_Nx_minusj = ps.ij_to_rc(Nx, j, Nx - 1, j, eq, 1);
            auto id_psi_Nxj =       ps.ij_to_rc(Nx, j, Nx,     j, eq, 1);
            A_nonzero.push_back(triplet(id_wNxj.first,          id_wNxj.second,           1.0));
            A_nonzero.push_back(triplet(id_psi_Nxj.first,       id_psi_Nxj.second,        two_div_hx_sq));
            A_nonzero.push_back(triplet(id_psi_Nx_minusj.first, id_psi_Nx_minusj.second, -two_div_hx_sq));
        };
    };

    // template  wp + apE * psi_E + app * psi_p + apW * psi_W + apS * psi_W + apN * psi_N = 0
    void build_system_laplace(const problem_params& ps, std::vector<triplet>& A_nonzero, const vector& sol_prev,
                              const vector& velocity_u, const vector& velocity_v, vector& rhs) {
        const size_type eq = 1; //here we are working with bot part of the matrix
        const size_type Nx = ps._Nx, Ny = ps._Ny; 
        const double div_tau = 1.0 / ps._tau, div_hx = 1.0 / ps._hx, div_hy = 1.0 / ps._hy;
        const double div_hx_sq = div_hx * div_hx, div_hy_sq = div_hy * div_hy;
        const double apE = div_hx_sq;
        const double app = -2.0 * (div_hx_sq + div_hy_sq);
        const double apW = div_hx_sq;
        const double apS = div_hy_sq;
        const double apN = div_hy_sq;

        // iterations over all internal grid nodes -> regular scheme template
        for (size_type i = 1; i < Nx; ++i)
            for (size_type j = 1; i < Ny; ++i) {
                //rhs
                size_type k = ps.ij_to_k(i, j, 1);
                rhs[k] = 0.0;

                //template
                auto id_wp = ps.ij_to_rc(i, j, i, j, eq, 0);
                A_nonzero.push_back(triplet(id_wp.first, id_wp.second, 1.0));

                auto id_psi_p = ps.ij_to_rc(i, j, i, j, eq, 1);
                A_nonzero.push_back(triplet(id_psi_p.first, id_psi_p.second, app));

                auto id_psi_E = ps.ij_to_rc(i, j, i + 1, j, eq, 1);
                A_nonzero.push_back(triplet(id_psi_E.first, id_psi_E.second, apE));

                auto id_psi_W = ps.ij_to_rc(i, j, i - 1, j, eq, 1);
                A_nonzero.push_back(triplet(id_psi_W.first, id_psi_W.second, apW));

                auto id_psi_S = ps.ij_to_rc(i, j, i, j + 1, eq, 1);
                A_nonzero.push_back(triplet(id_psi_S.first, id_psi_S.second, apS));

                auto id_psi_N = ps.ij_to_rc(i, j, i, j - 1, eq, 1);
                A_nonzero.push_back(triplet(id_psi_N.first, id_psi_N.second, apN));
            };

        // iterations over all boundary grid nodes -> boundary condition psi = 0
        const size_type count_x = ps._count_x, count_y = ps._count_y;
        for (size_type i = 0; i < count_x; ++i)     {   // BC bot
            auto id_psi_i0 = ps.ij_to_rc(i, 0, i, 0, eq, 1);
            A_nonzero.push_back(triplet(id_psi_i0.first, id_psi_i0.second, 1.0));
        };

        for (size_type j = 1; j < count_y - 1; ++j) {   // BC left (without corners) 
            auto id_psi_0j = ps.ij_to_rc(0, j, 0, j, eq, 1);
            A_nonzero.push_back(triplet(id_psi_0j.first, id_psi_0j.second, 1.0));
        };

        for (size_type i = 0; i < count_x; ++i)     {   // BC top
            auto id_psi_iNy = ps.ij_to_rc(i, Ny, i, Ny, eq, 1);
            A_nonzero.push_back(triplet(id_psi_iNy.first, id_psi_iNy.second, 1.0));
        };

        for (size_type j = 1; j < count_y - 1; ++j) {   // BC right (without corners)
            auto id_psi_Nxj = ps.ij_to_rc(Nx, j, Nx, j, eq, 1);
            A_nonzero.push_back(triplet(id_psi_Nxj.first, id_psi_Nxj.second, 1.0));
        };

    };

    void solve(const problem_params& ps, vector& sol_prev, vector& sol_curr, vector& velocity_u, vector& velocity_v, vector& rhs) {
        calc_velocity_on_vertices(ps, sol_prev, velocity_u, velocity_v);
        std::vector<triplet> A_nonzero;   
        A_nonzero.reserve(ps._size);                          
        build_system_transfer(ps, A_nonzero, sol_prev, velocity_u, velocity_v, rhs);
        build_system_laplace (ps, A_nonzero, sol_prev, velocity_u, velocity_v, rhs);
        sparse_LU(A_nonzero, sol_curr, rhs, false);
        sol_curr.swap(sol_prev);
    };
};