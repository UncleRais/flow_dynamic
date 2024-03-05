#include "library/vorcity_transfer.h"
#include "library/thermal_conductivity.h"
#include "library/proto_utils.hpp" /// ?


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

    std::ofstream log("log_integrals.txt");

    utl::progressbar::Percentage bar;
    bar.start();

    utl::timer::start();

    // Iterate over time
    for (uint step = 0; step < time_steps; ++step) {
        vorcity_transfer::solve(ps, sol_curr, sol_prev, velocity_u, velocity_v, rhs);

        // Export integrals related to conservation laws
        constexpr std::streamsize width = 16;
        log
            << std::setw(width) << vorcity_transfer::get_integral_omega(ps, sol_prev)
            << std::setw(width) << vorcity_transfer::get_integral_omega2(ps, sol_prev)
            << std::setw(width) << vorcity_transfer::get_integral_epsilon(ps, velocity_u, velocity_v)
            << std::setw(width) << vorcity_transfer::get_integral_dw_dt_psi(ps, sol_prev, sol_curr)
            << std::endl;

        bar.set_progress(static_cast<double>(step) / time_steps);
    }

    bar.finish();

    std::cout << "Completed in " << utl::timer::elapsed_string_fullform() << std::endl;

    utl::timer::start();

    ps.export_results(sol_prev, velocity_u, velocity_v);

    std::cout << "Exported in " << utl::timer::elapsed_string_fullform() << std::endl;

    return sol_curr;
}

int main(int argc, char** argv) {
    constexpr uint Nx = 100;
    constexpr uint Ny = 100;
    constexpr uint steps = 20;
    constexpr T L = 1.0;
    constexpr T H = 1.0;
    constexpr T time = 1.0;
    constexpr T nu = 1.0e6; // water: 1.787e6
    constexpr std::array<T, 4> boundary_velocities = { 0.0, 0.0, -1.0, -0.0 };
    constexpr auto filename = "result";
    constexpr auto format = SaveFormat::RAW;

    const problem_params params(Nx, Ny, steps, L, H, time, nu, boundary_velocities, filename, format);

    std::cout << utl::timer::datetime_string() << std::endl;

    const Vector solution = solve(params);

    return 0;
}