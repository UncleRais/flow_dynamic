#include "vorcity_transfer.h"
#include "thermal_conductivity.h"

// mat_coeffs: matrix nonzero elements
// sol: first half of elements - omega, second - psi
vector solve(const problem_params& ps) {
    size_type n = ps._size, time_steps = ps._steps;  
    vector rhs(n), sol_prev(n), sol_curr(n), velocity_u(ps._count_k), velocity_v(ps._count_k);

    // if needed, here we can work with initial T, psi and omega distributions
    // by default psi, omega = 0
    // T0 = delta_function(L / 2, H / 2) in order to obtain disturbance -> vorcity

    //time iterations
    for (size_type k = 0; k < time_steps; ++k) {
        vorcity_transfer::solve(ps, sol_curr, sol_prev, velocity_u, velocity_v, rhs);

    };
    
    return sol_curr;
};

void write_to_vtu(std::string filename, const vector& sol){

};


int main(int argc, char** argv) {
    // Nx = 2, Ny = 2, steps = 100, L = 1.0, H = 1.0, time = 1.0, nu = 1.0
    const problem_params params(2, 2, 100, 1.0, 1.0, 1.0, 1.0);
    vector sol = solve(params);
    std::cout << sol << std::endl;
    write_to_vtu("sol.vtu", sol);
    return 0;
};