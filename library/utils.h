#pragma once
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <optional>

#include <omp.h>

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "vtu11/vtu11.hpp"



using T         = double;
using size_type = std::size_t;
using Matrix    = Eigen::SparseMatrix<T>;
using Vector    = Eigen::VectorXd;
using Triplet   = Eigen::Triplet<T>;

struct MatrixIndex {
    size_type row;
    size_type col;
};

enum class Equation : size_type {
    FIRST = 0,
    SECOND = 1
};

enum class Variable : size_type {
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
    size_type _Nx;
    size_type _Ny;
    size_type _steps;

    T _L;
    T _H;
    T _time;

    T _nu;
    std::array<T, 4> _boundary_velocities;

    std::string _filename;
    std::string _folder = "results";
    SaveFormat _format;
        
    // 'Derived' parameters
    size_type _count_x;
    size_type _count_y;
    size_type _count_k;

    size_type _size;

    T  _hx;
    T  _hy;
    T  _tau;

    problem_params(
        size_type Nx, size_type Ny, size_type steps,
        T L, T H, T time,
        T nu, 
        std::array<T, 4> boundary_velocities,
        std::string_view filename,
        SaveFormat format
    ): 
        // 'Defining' parameters
        _Nx(Nx), _Ny(Ny), _steps(steps),
        _L(L), _H(H), _time(time),
        _nu(nu),
        _boundary_velocities{boundary_velocities},
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
        set_output_directory();
        _filename = _folder + "/" + _filename;
    }

    ~problem_params() = default;

    // mesh indexation convertation to solution vector indexation
    inline size_type ij_to_k(size_type i, size_type j) const {
        return i + j * _count_x;
    }

    // mesh indexation convertation to matrix indexation
    MatrixIndex ij_to_rc(
        size_type template_i,   size_type template_j,
        size_type term_i,       size_type term_j, 
        Equation equation,      Variable variable
    ) const {
        // Row corresponds to the 'k' of main template vertex
        const size_type row = ij_to_k(template_i, template_j) + static_cast<size_type>(equation) * _count_k;

        // Col corresponds to the 'k' of the term we want to add
        const size_type col = ij_to_k(term_i, term_j) + static_cast<size_type>(variable) * _count_k;

        return MatrixIndex{ row, col };
    }

    size_type ij_to_row(
        size_type template_i, size_type template_j,
        Equation equation
    ) const {
        // Row corresponds to the 'k' of main template vertex
        const size_type row = ij_to_k(template_i, template_j) + static_cast<size_type>(equation) * _count_k;

        return row;
    }

    std::string add_folder(std::string filename) const {
        return _folder + "/" + filename;
    }

    std::string str_i (size_type i) const {
        std::string str_i = std::to_string(i);
        while (str_i.size() < 6) str_i = "0" + str_i;
        return str_i;
    };
   
    void set_output_directory() const {
        std::filesystem::path path = "./" + _folder;
        if(!std::filesystem::is_directory(path))
            std::filesystem::create_directory(path);
        assert(std::filesystem::is_directory(path));
    }

    void export_results(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v, 
                        std::optional<std::string> filename) const {
        if (!filename) *filename = _filename; 
        if (_format == SaveFormat::VTU) {
            this->export_results_as_vtu(sol_prev, velocity_u, velocity_v, filename);
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

        for (size_type i = 0; i <= _Nx; ++i)
            for (size_type j = 0; j <= _Ny; ++j)
                file
                    << std::setw(width) << i * _hx
                    << std::setw(width) << j * _hy
                    << std::endl;

        file.close();

        // Omega
        file.open(_filename + "{omega}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (size_type i = 0; i <= _Nx; ++i)
            for (size_type j = 0; j <= _Ny; ++j)
                file
                    << std::setw(width) << sol_prev[ij_to_row(i, j, Equation::FIRST)]
                    << std::endl;

        file.close();

        // Psi
        file.open(_filename + "{psi}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (size_type i = 0; i <= _Nx; ++i)
            for (size_type j = 0; j <= _Ny; ++j)
                file
                << std::setw(width) << sol_prev[ij_to_row(i, j, Equation::SECOND)]
                << std::endl;

        file.close();

        // U + V
        file.open(_filename + "{uv}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (size_type i = 0; i <= _Nx; ++i)
            for (size_type j = 0; j <= _Ny; ++j)
                file
                    << std::setw(width) << velocity_u[ij_to_k(i, j)]
                    << std::setw(width) << velocity_v[ij_to_k(i, j)]
                    << std::endl;

        file.close();
    }

    void export_results_as_vtu(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v, 
                               std::optional<std::string> filename) const {

        std::vector<T> points;
        points.reserve(3 * _size);

        // Vertices
        for (size_type j = 0; j < _count_y; ++j) {
            for (size_type i = 0; i < _count_x; ++i) {
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

        size_type offset = 0;
        size_type type = 8;
        for (size_type j = 0; j < _Ny; ++j) {
            for (size_type i = 0; i < _Nx; ++i) {
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
        for (size_type k = 0; k < _count_k; ++k) {
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
        if (filename.has_value()) { 
            vtu11::writeVtu( *filename + ".vtu", mesh, dataSetInfo, { omega, psi, uv }, "Ascii" );
        } else {
            vtu11::writeVtu( _filename + ".vtu", mesh, dataSetInfo, { omega, psi, uv }, "Ascii" );
        }
    }

    void export_vtu_series(std::string filename) const {
        std::ofstream log(add_folder(filename) + ".vtu.series");
        const std::string tab("  ");
        const std::string _2tab(tab + tab);

        log << "{" << std::endl;
        log << tab + "\"file-series-version\" : \"1.0\"," << std::endl; 
        log << tab + "  \"files\" : [" << std::endl;
        for (size_type k; k < _steps; ++k)
            log << _2tab + "{ \"name\" : \"result_" << str_i(k) << ".vtu\", \"time\" : " << k * _tau << " }," << std::endl;
        log << _2tab + "{ \"name\" : \"result_" << str_i(_steps) << ".vtu\", \"time\" : " << _steps * _tau << " }" << std::endl;
        log << "  ]" << std::endl;
        log << "}" << std::endl;
    }
};

std::string system_as_string(const Matrix &A, const Vector &rhs) {
    const size_type size = A.rows();

    assert(size == A.rows());
    assert(size == A.cols());
    assert(size == rhs.rows());

    std::stringstream ss;

    ss << "SLAE [" << size << " x "<<  size << "]:" << std::endl;

    for (size_type i = 0; i < size; ++i) {
        ss << "[ ";

        for (size_type j = 0; j < size; ++j) {

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