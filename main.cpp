#include "library/test_methods.h"
#include "library/vorcity_transfer.h"
#include "library/thermal_conductivity.h"
#include "library/proto_utils.hpp" 

Vector test_128() {
    constexpr size_type Nx = 128;
    constexpr size_type Ny = 128;
    constexpr size_type steps = 100;
    constexpr T L = 1.0;
    constexpr T H = 1.0;
    constexpr T time = 100.0;
    constexpr T Re = 1000.0;
    constexpr std::array<T, 4> boundary_velocities = { 0.0, 0.0, 1.0, 0.0 };
    constexpr auto filename = "result";
    constexpr auto format = SaveFormat::RAW;
    constexpr Method method = DecompositionMethod::SPARSE_LU;
    //constexpr Method method = IterativeMethodData{ IterativeMethod::BI_CGSTAB, 1e-12, 1000 };

    const problem_params ps(Nx, Ny, steps, L, H, time, 1. / Re, boundary_velocities, method, filename, format);

    const size_type size = ps._size;
    const size_type time_steps = ps._steps;

    Vector rhs = Vector::Zero(size);
    Vector sol_prev = Vector::Zero(size);
    Vector sol_curr = Vector::Zero(size);
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

    std::ofstream log(ps.get_relative_path("log_integrals.txt"));

    std::cout
        << "Running 'test_128()'..." << std::endl
        << utl::timer::datetime_string() << std::endl;

    utl::progressbar::Percentage bar;
    bar.start();
    utl::timer::start();

    // Iterate over time
    for (size_type step = 0; step < time_steps; ++step) {
        vorcity_transfer::solve(ps, sol_curr, sol_prev, velocity_u, velocity_v, rhs);

        // Export integrals related to conservation laws
        constexpr std::streamsize width = 16;
        log
            << std::setw(width) << test::get_integral_omega(ps, sol_prev)
            << std::setw(width) << test::get_integral_omega2(ps, sol_prev)
            << std::setw(width) << test::get_integral_epsilon(ps, velocity_u, velocity_v)
            << std::setw(width) << test::get_integral_dw_dt_psi(ps, sol_prev, sol_curr)
            << std::endl;

        bar.set_progress(static_cast<double>(step) / time_steps);
    }

    bar.finish();
    std::cout << "Completed in " << utl::timer::elapsed_string_fullform() << std::endl;
    
    ps.export_results(sol_prev, velocity_u, velocity_v);

    // Test against paper results "High-Re Solutions for Incompressible Flow by U. GHIA"
    constexpr std::array test_i { 0, 8, 30, 64, 103, 124, 128 };
    constexpr std::array test_v { 0.00000, 0.27485, 0.32235, 0.02526, -0.31966, -0.21388, 0.00000 };

    for (size_type testcase = 0; testcase < test_i.size(); ++testcase) {
        const auto i = test_i[testcase];
        const auto j = 64;
        const auto k = ps.ij_to_k(i, j);

        const T v_test   = test_v[testcase];
        const T v_actual = velocity_v(k);

        const T error = std::abs(v_test - v_actual);

        std::cout << "error for v[" << i << ", " << j << "]: " << error << std::endl;
    }

    return sol_curr;
}

Vector solve(const problem_params &ps, bool file_logging = true, bool tests = true) {
    const size_type size        = ps._size;
    const size_type time_steps  = ps._steps;
    const SaveFormat format     = ps._format;

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

    std::ofstream log(ps.get_relative_path("log_integrals.txt"));

    std::cout
        << "Running 'solve()'..." << std::endl
        << utl::timer::datetime_string() << std::endl;

    utl::progressbar::Percentage bar;
    bar.start();

    utl::timer::start();
    const auto filename = ps.get_relative_path("result_");

    if (file_logging) ps.export_results(sol_prev, velocity_u, velocity_v, filename + ps.str_i(0));

    // Iterate over time
    for (size_type step = 0; step < time_steps; ++step) {
        vorcity_transfer::solve(ps, sol_curr, sol_prev, velocity_u, velocity_v, rhs);

        if (file_logging && format == SaveFormat::VTU)
            ps.export_results(sol_prev, velocity_u, velocity_v, filename + ps.str_i(step + 1));
            // we only export intermediate time layers if it's VTU series

        // Export integrals related to conservation laws
        constexpr std::streamsize width = 16;
        if (tests) log
            << std::setw(width) << test::get_integral_omega(ps, sol_prev)
            << std::setw(width) << test::get_integral_omega2(ps, sol_prev)
            << std::setw(width) << test::get_integral_epsilon(ps, velocity_u, velocity_v)
            << std::setw(width) << test::get_integral_dw_dt_psi(ps, sol_prev, sol_curr)
            << std::endl;

        bar.set_progress(static_cast<double>(step) / time_steps);
    }
    
    bar.finish();
    std::cout << "Completed in " << utl::timer::elapsed_string_fullform() << std::endl;

    if (file_logging && format == SaveFormat::VTU) ps.export_vtu_series_metadata("timesteps");
    else                                           ps.export_results(sol_prev, velocity_u, velocity_v);
    
    return sol_curr;
}

int main(int argc, char** argv) {
    constexpr size_type Nx = 200;
    constexpr size_type Ny = 200;
    constexpr size_type steps = 200;
    constexpr T L = 1.0;
    constexpr T H = 1.0;
    constexpr T time = 100.0;
    constexpr T nu = 1e-3; // water: 1.787e-6
    constexpr std::array<T, 4> boundary_velocities = { 0.0, 0.0, -1.0, -0.0 };
    constexpr auto filename = "result";
    constexpr auto format = SaveFormat::RAW;
    constexpr Method method = DecompositionMethod::SPARSE_LU;
    //constexpr Method method = IterativeMethodData{ IterativeMethod::BI_CGSTAB, 1e-12, 1000 };

    const problem_params params(Nx, Ny, steps, L, H, time, nu, boundary_velocities, method, filename, format);

    //const Vector solution = solve(params, true, true);
    test_128();

    return 0;
}