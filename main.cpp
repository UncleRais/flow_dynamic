#include "library/vorcity_transfer.h"
#include "library/thermal_conductivity.h"


Vector solve(const problem_params &ps) {
    const uint size        = ps._size;
    const uint time_steps = ps._steps;

    Vector rhs        = Vector::Zero(size);
    Vector sol_prev   = Vector::Zero(size);
    Vector sol_curr   = Vector::Zero(size);
    Vector velocity_u = Vector::Zero(ps._count_k);
    Vector velocity_v = Vector::Zero(ps._count_k);
        // sol[ 0      ... size/2 ] corresponds to 'omega'
        // sol[ size/2 ... size   ] corresponds to 'psi'

    // Initial conditions:
    // { psi   |t=0   = 0
    // { omega |t=0   = 0
    // { T     |t=0   = dirac_delta(L/2, H/2)
    // The choise of initial 'T' introduces disturbance to the system => vorcity
    vorcity_transfer::set_velocity_on_boundary(ps, velocity_u, velocity_v);

    // Iterate over time
    for (uint step = 0; step < time_steps; ++step) {
        vorcity_transfer::solve(ps, sol_curr, sol_prev, velocity_u, velocity_v, rhs);
    }

    ps.export_results(sol_prev, velocity_u, velocity_v);

    return sol_curr;
}

int main(int argc, char** argv) {
    constexpr uint Nx = 3;
    constexpr uint Ny = 2;
    constexpr uint steps = 10;
    constexpr T L = 1.0;
    constexpr T H = 1.0;
    constexpr T time = 1.0;
    constexpr T nu = 1.0; // water: 1.787e6
    constexpr std::array<T, 4> boundary_velocities = { 1.0, 0.0, 0.0, 0.0 };
    constexpr auto filename = "result";
    constexpr auto format = SaveFormat::RAW;

    const problem_params params(Nx, Ny, steps, L, H, time, nu, boundary_velocities, filename, format);

    const Vector solution = solve(params);

    return 0;
}