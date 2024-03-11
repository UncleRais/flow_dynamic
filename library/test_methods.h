#pragma once

#include <algorithm>

#include "utils.h"
#include "slae_solver.h"


namespace test {

    inline T integrate_value_over_G(const problem_params &ps, const Vector &values, std::function<T(T)> func) {
        const size_type Nx = ps._Nx;
        const size_type Ny = ps._Ny;

        const size_type size = ps._size;
        const T area = ps._hx * ps._hy;

        // Compure epsilon = 0.5 * integral[u^2 + v^2, G]
        T sum = 0.0;

        for (size_type i = 0; i <= Nx - 1; ++i) {
            for (size_type j = 0; j <= Ny - 1; ++j) {
                const auto idx_bl = ps.ij_to_k(i,     j    );
                const auto idx_br = ps.ij_to_k(i + 1, j    );
                const auto idx_tr = ps.ij_to_k(i + 1, j + 1);
                const auto idx_tl = ps.ij_to_k(i,     j + 1);

                const T val_bl = values[idx_bl];
                const T val_br = values[idx_br];
                const T val_tr = values[idx_tr];
                const T val_tl = values[idx_tl];

                sum += (func(val_bl) + func(val_br) + func(val_tr) + func(val_tl)) * area * 0.25;
            }
        }

        return sum;
    }


    inline T get_integral_omega(const problem_params &ps, const Vector &sol) {
        const size_type half_size = ps._size / 2;

        // Extract omega from solution
        Vector omega(half_size);
        std::copy(sol.begin(), sol.begin() + half_size, omega.begin());

        const auto identity = [](T x) -> T { return x; };

        // Compute sum = integral[omega, G]
        const T sum = integrate_value_over_G(ps, omega, identity);

        return sum;
    }


    inline T get_integral_omega2(const problem_params &ps, const Vector &sol) {
        const size_type half_size = ps._size / 2;

        // Extract omega from solution
        Vector omega(half_size);
        std::copy(sol.begin(), sol.begin() + half_size, omega.begin());

        const auto sqr = [](T x) -> T { return x * x; };

        // Compute sum = integral[omega, G]
        const T sum = integrate_value_over_G(ps, omega, sqr);

        return sum;
    }


    inline T get_integral_dw_dt_psi(const problem_params &ps, const Vector &sol_prev, const Vector &sol_next) {
        const size_type half_size = ps._size / 2;

        // Extract psi, omega from solution
        Vector psi_prev(half_size);
        std::copy(sol_prev.begin() + half_size, sol_prev.end(), psi_prev.begin());

        Vector psi_next(half_size);
        std::copy(sol_prev.begin() + half_size, sol_prev.end(), psi_next.begin());

        Vector omega_prev(half_size);
        std::copy(sol_prev.begin(), sol_prev.begin() + half_size, omega_prev.begin());

        Vector omega_next(half_size);
        std::copy(sol_next.begin(), sol_next.begin() + half_size, omega_next.begin());

        Vector dw_dt_psi = Vector::Zero(half_size);
        for (size_type i = 0; i < half_size; ++i)
            dw_dt_psi[i] += (omega_next[i] - omega_prev[i]) / ps._tau * 0.5 * (psi_prev[i] + psi_next[i]);

        const auto identity = [](T x) -> T { return x; };

        // Compute sum = integral[omega, G]
        const T I = integrate_value_over_G(ps, dw_dt_psi, identity);

        return I;
    }


    inline T get_integral_epsilon(const problem_params &ps, const Vector &velocity_u, const Vector &velocity_v) {

        // Compure epsilon = 0.5 * integral[u^2 + v^2, G]
        const auto sqr = [](T x) -> T { return x * x; };

        const T epsilon = 0.5 * (integrate_value_over_G(ps, velocity_u, sqr) + integrate_value_over_G(ps, velocity_v, sqr));

        return epsilon;
    }

    inline void build_test_problem(std::vector<Triplet> &coefficients, Eigen::VectorXd &b, size_type n){
    for (size_type i = 0; i < n; ++i) {
        b[i] = 0.0001;
        for (size_type j = 0; j < n; ++j) 
            if(i == j) {
                const T val = (j + 1) / 1000.0;
                coefficients.push_back(Triplet(i, j, val));
            };
    };
};
}