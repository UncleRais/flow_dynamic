#include "library/test_methods.h"
#include "library/vorcity_transfer.h"
#include "library/thermal_conductivity.h"
#include "library/proto_utils.hpp" 

Vector test_128() {

    // REFERENCE solution
    // Test against paper results "High-Re Solutions for Incompressible Flow by U. GHIA"
    Vector test_i(17);
           test_i <<   0,   8,   9,   10, 
                       12,  20,  29,  30, 
                       64,  103, 110, 116, 
                       121, 122, 123, 124, 128 ;
    Vector test_v(17);
           test_v << 0.00000,  0.27485,  0.29012,  0.30353, 
                     0.32627,  0.37095,  0.33075,  0.32235,
                     0.02526, -0.31966, -0.42665, -0.51550,
                    -0.39188, -0.33714, -0.27669, -0.21388,  0.00000;

    auto compare_to_reference = [&](const problem_params& ps, const Vector& velocity_v, bool verbose = false) {
        // Choose necessary values 
        Vector actual_v = Vector::Zero(test_v.size());
        std::vector<std::pair<size_type, size_type>> actual_ij(test_i.size(), std::make_pair(0, 0));
        for (size_type testcase = 0; testcase < test_i.size(); ++testcase) { 
            auto& ij = actual_ij[testcase];
            ij.first = test_i[testcase];
            ij.second = 64;
            const auto k = ps.ij_to_k(ij);
            actual_v[testcase] = velocity_v(k);
        }

        // Calculate point-to-point relative error
        if (verbose) relative_error(test_v, actual_v, actual_ij);
        return (actual_v - test_v).norm() / (test_v.norm() + 1e-16);
    };

    constexpr size_type Nx = 128;
    constexpr size_type Ny = 128;
    constexpr size_type steps = 100;
    constexpr T L = 1.0;
    constexpr T H = 1.0;
    constexpr T time = 100.0;
    constexpr T Re = 1000.0;
    constexpr T kappa = 0.0;
    constexpr T lambda = 1e-3;
    constexpr T Ra = 0.0;
    constexpr std::array<T, 4> boundary_velocities = { 0.0, 0.0, 1.0, 0.0 };
    constexpr std::array<std::pair<ThermalBoundaryType, T>, 4> boundary_temperature_conditions 
                                                        = { std::make_pair(ThermalBoundaryType::TEMPERATURE, 0.0), 
                                                            std::make_pair(ThermalBoundaryType::TEMPERATURE, 0.0), 
                                                            std::make_pair(ThermalBoundaryType::TEMPERATURE, 0.0),  
                                                            std::make_pair(ThermalBoundaryType::TEMPERATURE, 0.0) };
    const std::function<T(T, T)> initial_temperature = [&](T x, T y) { return T(0.0);};

    constexpr auto filename = "result";
    constexpr auto format = SaveFormat::RAW;
    constexpr Method method = DecompositionMethod::SPARSE_LU;
    //constexpr Method method = IterativeMethodData{ IterativeMethod::BI_CGSTAB, 1e-12, 1000 };

    const problem_params ps(Nx, Ny, steps, L, H, time, 1. / Re, kappa, lambda, Ra, 
                            boundary_velocities,boundary_temperature_conditions, initial_temperature, method, filename, format);

    const size_type size = ps._size;
    const size_type time_steps = ps._steps;

    Vector rhs = Vector::Zero(size);
    Vector sol_prev = Vector::Zero(size);
    Vector sol_curr = Vector::Zero(size);
    Vector temperature = Vector::Zero(ps._count_k);
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
    constexpr std::streamsize width = 20;

    log
            << std::setw(width) << "omega"
            << std::setw(width) << "omega^2"
            << std::setw(width) << "epsilon"
            << std::setw(width) << "dw/dt*psi"
            << std::setw(width) << "div{u, v, 0}"
            << std::setw(width) << "relative_err"
            << std::endl;

    // Iterate over time
    for (size_type step = 0; step < time_steps; ++step) {
        vorcity_transfer::solve(ps, sol_curr, sol_prev, velocity_u, velocity_v, temperature, rhs);

        // Export integrals related to conservation laws
        log
            << std::setw(width) << test::get_integral_omega(ps, sol_prev)
            << std::setw(width) << test::get_integral_omega2(ps, sol_prev)
            << std::setw(width) << test::get_integral_epsilon(ps, velocity_u, velocity_v)
            << std::setw(width) << test::get_integral_dw_dt_psi(ps, sol_prev, sol_curr)
            << std::setw(width) << test::get_div_uv(ps, velocity_u, velocity_v)                      
            << std::setw(width) << compare_to_reference(ps, velocity_v, false)
            << std::endl;

        bar.set_progress(static_cast<double>(step) / time_steps);
    }

    bar.finish();
    std::cout << "Completed in " << utl::timer::elapsed_string_fullform() << std::endl;
    
    ps.export_results(sol_prev, temperature, velocity_u, velocity_v);

    auto error = compare_to_reference(ps, velocity_v, true);
    std::cout << "Relative vector error = " << error << std::endl;

    return sol_curr;
}

Vector solve(const problem_params &ps, bool file_logging = true, bool tests = true) {
    const size_type size        = ps._size;
    const size_type time_steps  = ps._steps;
    const SaveFormat format     = ps._format;

    Vector rhs_vorcity = Vector::Zero(size);
    Vector sol_prev    = Vector::Zero(size);
    Vector sol_curr    = Vector::Zero(size);
    Vector rhs_temp         = Vector::Zero(ps._count_k);
    Vector temperature_prev = Vector::Zero(ps._count_k);
    Vector temperature_curr = Vector::Zero(ps._count_k);
    Vector velocity_u = Vector::Zero(ps._count_k);
    Vector velocity_v = Vector::Zero(ps._count_k);
    // sol[ 0      ... size/2 ] corresponds to 'omega'
    // sol[ size/2 ... size   ] corresponds to 'psi'

    // Initial conditions:
    // { psi   |t=0   = 0
    // { omega |t=0   = 0
    // { T     |t=0   = <Linear distribution>
    // The choise of initial 'T' introduces disturbance to the system => vorcity
    vorcity_transfer::set_velocity_on_boundary(ps, velocity_u, velocity_v);
    thermal_conductivity::set_initial_temperature(ps, temperature_prev);

    std::ofstream log(ps.get_relative_path("log_integrals.txt"));

    constexpr std::streamsize width = 20;
    log
            << std::setw(width) << "omega"
            << std::setw(width) << "omega^2"
            << std::setw(width) << "epsilon"
            << std::setw(width) << "dw/dt*psi"
            << std::setw(width) << "div{u, v, 0}"
            << std::endl;

    std::cout
        << "Running 'solve()'..." << std::endl
        << utl::timer::datetime_string() << std::endl;

    utl::progressbar::Percentage bar;
    bar.start();

    utl::timer::start();
    const auto filename = ps.get_relative_path("result_");

    if (file_logging) ps.export_results(sol_prev, temperature_prev, velocity_u, velocity_v, filename + ps.str_i(0));

    // Iterate over time
    for (size_type step = 0; step < time_steps; ++step) {
        vorcity_transfer::solve    (ps, sol_prev,         sol_curr,         velocity_u, velocity_v, temperature_prev, rhs_vorcity);
        thermal_conductivity::solve(ps, temperature_prev, temperature_curr, velocity_u, velocity_v, rhs_temp);

        if (file_logging && format == SaveFormat::VTU)
            ps.export_results(sol_prev, temperature_prev, velocity_u, velocity_v, filename + ps.str_i(step + 1));
            // we only export intermediate time layers if it's VTU series

        // Export integrals related to conservation laws
        if (tests) log
            << std::setw(width) << test::get_integral_omega(ps, sol_prev)
            << std::setw(width) << test::get_integral_omega2(ps, sol_prev)
            << std::setw(width) << test::get_integral_epsilon(ps, velocity_u, velocity_v)
            << std::setw(width) << test::get_integral_dw_dt_psi(ps, sol_prev, sol_curr)
            << std::setw(width) << test::get_div_uv(ps, velocity_u, velocity_v)
            << std::endl;

        bar.set_progress(static_cast<double>(step) / time_steps);
    }
    
    bar.finish();
    std::cout << "Completed in " << utl::timer::elapsed_string_fullform() << std::endl;

    if (file_logging && format == SaveFormat::VTU) ps.export_vtu_series_metadata("timesteps");
    else                                           ps.export_results(sol_prev, temperature_prev, velocity_u, velocity_v);
    
    return sol_curr;
}

int main(int argc, char** argv) {
    constexpr size_type Nx = 26;
    constexpr size_type Ny = 26;
    constexpr size_type steps = 2000;
    constexpr T L = 1.0;
    constexpr T H = 1.0;
    constexpr T time = 200.0;
    constexpr T nu = 1; // water: 1.787e-6
    constexpr T kappa = 1;
    constexpr T lambda = 1;
    constexpr T Ra = 30080;
    constexpr std::array<T, 4> boundary_velocities = { 0.0, 0.0, 0.0, 0.0 };
    constexpr std::array<std::pair<ThermalBoundaryType, T>, 4> boundary_temperature_conditions 
                                                        = { std::make_pair(ThermalBoundaryType::TEMPERATURE, 0.0), 
                                                            std::make_pair(ThermalBoundaryType::FLUX, 0.0), 
                                                            std::make_pair(ThermalBoundaryType::TEMPERATURE, 1.0),  
                                                            std::make_pair(ThermalBoundaryType::FLUX, 0.0) };
    constexpr T var = 1 / L * 2 * M_PI;
    const std::function<T(T, T)> initial_temperature = [&](T x, T y) { 
        if (std::abs(y - H * 0.33333333) < 0.01) 
            return y + 0.01 * std::sin(x * var);
        return y;
    };

    constexpr auto filename = "result";
    constexpr auto format = SaveFormat::VTU;
    constexpr Method method = DecompositionMethod::SPARSE_LU;
    //constexpr Method method = IterativeMethodData{ IterativeMethod::BI_CGSTAB, 1e-12, 1000 };

    const problem_params params(Nx, Ny, steps, L, H, time, nu, kappa, lambda, Ra,
                                boundary_velocities, boundary_temperature_conditions, initial_temperature, method, filename, format);

    const Vector solution = solve(params, true, true);
    //test_128();

    return 0;
}