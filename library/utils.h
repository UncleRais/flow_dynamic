#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

#include <omp.h>

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "vtu11/vtu11.hpp"



using T         = double;
using uint      = size_t;
using Matrix    = Eigen::SparseMatrix<T>;
using Vector    = Eigen::VectorXd;
using Triplet   = Eigen::Triplet<T>;

struct MatrixIndex {
    uint row;
    uint col;
};

enum class Equation : uint {
    FIRST = 0,
    SECOND = 1
};

enum class Variable : uint {
    OMEGA = 0,
    PSI = 1
};

enum class SaveFormat {
    RAW,
    VTU
};

inline void exit_with_error(std::string_view message) {
    std::cerr << "ERROR: " << message << std::endl;
    exit(1);
}

// this struct contains all information about problem:
// 1) physical characteristics 2) computational domain sizes
// 3) mesh 4) boundary condition 5) initial distributions
 
struct problem_params {
    // 'Defining' parameters
    uint _Nx;
    uint _Ny;
    uint _steps;

    T _L;
    T _H;
    T _time;

    T _nu;
    std::array<T, 4> _boundary_velocities;

    std::string _filename;
    SaveFormat _format;
        
    // 'Derived' parameters
    uint _count_x;
    uint _count_y;
    uint _count_k;

    uint _size;

    T  _hx;
    T  _hy;
    T  _tau;

    problem_params(
        uint Nx, uint Ny, uint steps,
        T L, T H, T time,
        T nu, 
        std::array<T, 4> BC_velocity,
        std::string_view filename,
        SaveFormat format
    ): 
        // 'Defining' parameters
        _Nx(Nx), _Ny(Ny), _steps(steps),
        _L(L), _H(H), _time(time),
        _nu(nu),
        _boundary_velocities{BC_velocity},
        _filename(filename),
        _format(format)
    {
        // 'Derived' parameters
        _hx = L / Nx; 
        _hy = H / Ny; 
        _tau = time / steps;
        _count_x = _Nx + 1;
        _count_y = _Ny + 1;
        _count_k = _count_x * _count_y;
        _size = 2 * _count_k;
    }

    ~problem_params() = default;

    // mesh indexation convertation to solution vector indexation
    inline uint ij_to_k(uint i, uint j) const {
        return i + j * _count_x;
    }

    // mesh indexation convertation to matrix indexation
    MatrixIndex ij_to_rc(
        uint template_i,   uint template_j,
        uint term_i,       uint term_j, 
        Equation equation, Variable variable
    ) const {
        // Row corresponds to the 'k' of main template vertex
        const uint row = ij_to_k(template_i, template_j) + static_cast<uint>(equation) * _count_k;

        // Col corresponds to the 'k' of the term we want to add
        const uint col = ij_to_k(term_i, term_j) + static_cast<uint>(variable) * _count_k;

        return MatrixIndex{ row, col };
    }

    uint ij_to_row(
        uint template_i, uint template_j,
        Equation equation
    ) const {
        // Row corresponds to the 'k' of main template vertex
        const uint row = ij_to_k(template_i, template_j) + static_cast<uint>(equation) * _count_k;

        return row;
    }

    void export_results(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v) const {
        if (_format == SaveFormat::VTU) {
            this->export_results_as_vtu(sol_prev, velocity_u, velocity_v);
        }
        else if (_format == SaveFormat::RAW) {
            this->export_results_as_raw(sol_prev, velocity_u, velocity_v);
        }
        else {
            exit_with_error("Unrecognized save format.");
        }
    }

    void export_results_as_raw(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v) const {
        constexpr auto extension = ".txt";
        constexpr auto error     = "Could not create file.";
        constexpr std::streamsize width = 20;

        std::ofstream file;

        // NOTE:
        // We must maintain the same order of iteration (i->j) so each value correspond to the same vertex
        // No point in saving all indexation when exporting raw data

        // Vertices
        file.open(_filename + "{vertices}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (uint i = 0; i <= _Nx; ++i)
            for (uint j = 0; j <= _Ny; ++j)
                file
                    << std::setw(width) << i * _hx
                    << std::setw(width) << j * _hy
                    << std::endl;

        file.close();

        // Omega
        file.open(_filename + "{omega}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (uint i = 0; i <= _Nx; ++i)
            for (uint j = 0; j <= _Ny; ++j)
                file
                    << std::setw(width) << sol_prev[ij_to_row(i, j, Equation::FIRST)]
                    << std::endl;

        file.close();

        // Psi
        file.open(_filename + "{psi}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (uint i = 0; i <= _Nx; ++i)
            for (uint j = 0; j <= _Ny; ++j)
                file
                << std::setw(width) << sol_prev[ij_to_row(i, j, Equation::SECOND)]
                << std::endl;

        file.close();

        // U + V
        file.open(_filename + "{uv}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (uint i = 0; i <= _Nx; ++i)
            for (uint j = 0; j <= _Ny; ++j)
                file
                    << std::setw(width) << velocity_u[ij_to_k(i, j)]
                    << std::setw(width) << velocity_v[ij_to_k(i, j)]
                    << std::endl;

        file.close();
    }

    void export_results_as_vtu(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v) const {

        std::vector<T> points;
        points.reserve(3 * _size);

        // Vertices
        for (uint j = 0; j < _count_y; ++j) {
            for (uint i = 0; i < _count_x; ++i) {
                points.push_back(i * _hx);
                points.push_back(j * _hy);
                points.push_back(T(0.0));
            }
        }

        // Mesh geometry
        std::vector<vtu11::VtkIndexType> connectivity;
        std::vector<vtu11::VtkIndexType> offsets;
        std::vector<vtu11::VtkCellType> types;
        connectivity.reserve(_count_k);
        offsets.reserve(_count_k);
        types.reserve(_count_k);

        uint offset = 0;
        uint type = 8;
        for (uint j = 0; j < _Ny; ++j) {
            for (uint i = 0; i < _Nx; ++i) {
                offset += 4;
                connectivity.push_back(i + j * _count_x);
                connectivity.push_back(i + 1 + j * _count_x);
                connectivity.push_back(i + (j + 1) * _count_x);
                connectivity.push_back(i + 1 + (j + 1) * _count_x);
                offsets.push_back(offset);
                types.push_back(static_cast<vtu11::VtkCellType>(type));
            }
        }

        // Create small proxy mesh type
        vtu11::Vtu11UnstructuredMesh mesh { points, connectivity, offsets, types };

        // Omega, Psi, U + V
        std::vector<T> omega;
        std::vector<T> psi;
        std::vector<T> uv;
        omega.resize(_count_k);
        psi.resize(_count_k);
        std::copy(sol_prev.begin(),            sol_prev.begin() + _count_k, omega.begin());
        std::copy(sol_prev.begin() + _count_k, sol_prev.end(),              psi.begin());
        for (uint k = 0; k < _count_k; ++k) {
            uv.push_back(velocity_u[k]);
            uv.push_back(velocity_v[k]);
            uv.push_back(T(0.0));
        }

        // Create tuples with (name, association, number of components) for each data set
        std::vector<vtu11::DataSetInfo> dataSetInfo
        {
            { "Omega",      vtu11::DataSetType::PointData, 1 },
            { "Psi",        vtu11::DataSetType::PointData, 1 },
            { "Velocity",   vtu11::DataSetType::PointData, 3 }
        };

        // Write data to .vtu file using Ascii format
        vtu11::writeVtu( _filename + ".vtu", mesh, dataSetInfo, { omega, psi, uv }, "Ascii" );
    }
};

std::string system_as_string(const Matrix &A, const Vector &rhs) {
    const uint size = A.rows();

    assert(size == A.rows());
    assert(size == A.cols());
    assert(size == rhs.rows());

    std::stringstream ss;

    ss << "SLAE [" << size << " x "<<  size << "]:" << std::endl;

    for (uint i = 0; i < size; ++i) {
        ss << "[ ";

        for (uint j = 0; j < size; ++j) {

            const T coef = A.coeff(i, j);
            constexpr std::streamsize width = 9;

            if (coef == 0)
                ss << std::setw(width) << "-" << " ";
            else
                ss << std::setw(width) << std::setprecision(3) << coef << " ";
        }

        ss << " ]  {" << std::setw(10) << std::setprecision(3) << rhs.coeff(i) << "}"<< std::endl;
    }

    return ss.str();
}

std::ostream& operator<<(std::ostream &cout, const Vector &vec) {
    
    std::cout << "{ ";
    for (const auto &el: vec) std::cout << el << ", ";
    std::cout << " }";

    return cout;
}

/// What is that?
void build_test_problem(std::vector<Triplet> &coefficients, Eigen::VectorXd &b, uint n){
    for (uint i = 0; i < n; ++i) {
        b[i] = 0.0001;
        for (uint j = 0; j < n; ++j) 
            if(i == j) {
                const T val = (j + 1) / 1000.0;
                coefficients.push_back(Triplet(i, j, val));
            };
    };
};