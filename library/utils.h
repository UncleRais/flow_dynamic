#pragma once
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <optional>
#include <variant>

#include <omp.h>

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "vtu11/vtu11.hpp"



using T         = double;
using size_type = int;
using Matrix    = Eigen::SparseMatrix<T/*, Eigen::RowMajor*/>; /// for some reason row-major storage causes crashed for large N
using Vector    = Eigen::VectorXd;
using Triplet   = Eigen::Triplet<T>;

struct MatrixIndex {
    size_type row;
    size_type col;
};

// Indexation-related enums
enum class Equation : size_type {
    FIRST = 0,
    SECOND = 1
};

enum class Variable : size_type {
    OMEGA = 0,
    PSI = 1
};

// SLAE enums
enum class DecompositionMethod {
    SPARSE_LU,
    SPARSE_QR,
    SIMPLICIAL_CHOLESKY
};

enum class IterativeMethod {
    LS_CONJUGATE_GRADIENT,
    BI_CGSTAB
};

struct IterativeMethodData {
    IterativeMethod                               method;
    std::optional<T>                              tolerance      = std::nullopt; // default: ~2.2e-16
    std::optional<Eigen::Index>                   max_iterations = std::nullopt; // default:  2 * matrix_size
};

using Method = std::variant<DecompositionMethod, IterativeMethodData>;

// Format enums
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

    Method _method;

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
        Method method,
        std::string_view filename,
        SaveFormat format
    ): 
        // 'Defining' parameters
        _Nx(Nx), _Ny(Ny), _steps(steps),
        _L(L), _H(H), _time(time),
        _nu(nu),
        _boundary_velocities(boundary_velocities),
        _method(method),
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
    inline MatrixIndex ij_to_rc(
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

    inline size_type ij_to_rhs_row(
        size_type template_i, size_type template_j,
        Equation equation
    ) const {
        // Row corresponds to the 'k' of main template vertex
        const size_type row = ij_to_k(template_i, template_j) + static_cast<size_type>(equation) * _count_k;

        return row;
    }

    inline size_type ij_to_sol_row(
        size_type template_i, size_type template_j,
        Variable variable
    ) const {
        // Row corresponds to the 'k' of main template vertex
        const size_type row = ij_to_k(template_i, template_j) + static_cast<size_type>(variable) * _count_k;

        return row;
    }

    inline std::string get_relative_path(std::string filename) const {
        return _folder + "/" + filename;
    }

    inline std::string str_i (size_type i) const {
        std::string str_i = std::to_string(i);
        while (str_i.size() < 6) str_i = "0" + str_i;
        return str_i;
    };
   
    inline void set_output_directory() const {
        std::filesystem::path path = "./" + _folder;
        if(!std::filesystem::is_directory(path))
            std::filesystem::create_directory(path);
        assert(std::filesystem::is_directory(path));
    }

    inline void export_results(
        const Vector &sol_prev,
        const Vector &velocity_u, const Vector &velocity_v, 
        std::optional<std::string> filename = std::nullopt
    ) const {
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

    inline void export_results_as_raw(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v) const {
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
                    << std::setw(width) << sol_prev[ij_to_sol_row(i, j, Variable::OMEGA)]
                    << std::endl;

        file.close();

        // Psi
        file.open(_filename + "{psi}" + extension);
        if (!file.is_open()) exit_with_error(error);

        for (size_type i = 0; i <= _Nx; ++i)
            for (size_type j = 0; j <= _Ny; ++j)
                file
                << std::setw(width) << sol_prev[ij_to_sol_row(i, j, Variable::PSI)]
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

    inline void export_results_as_vtu(const Vector &sol_prev, const Vector &velocity_u, const Vector &velocity_v,
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

    inline void export_vtu_series_metadata(std::string filename) const {
        std::ofstream log(get_relative_path(filename) + ".vtu.series");
        const std::string tab("  ");
        const std::string _2tab(tab + tab);

        log << "{" << std::endl;
        log << tab + "\"file-series-version\" : \"1.0\"," << std::endl; 
        log << tab + "  \"files\" : [" << std::endl;
        for (size_type k = 0; k < _steps; ++k)
            log << _2tab + "{ \"name\" : \"result_" << str_i(k) << ".vtu\", \"time\" : " << k * _tau << " }," << std::endl;
        log << _2tab + "{ \"name\" : \"result_" << str_i(_steps) << ".vtu\", \"time\" : " << _steps * _tau << " }" << std::endl;
        log << "  ]" << std::endl;
        log << "}" << std::endl;
    }
};

inline std::string system_as_string(const Matrix &A, const Vector &rhs) {
    const auto size = A.rows();

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

inline std::ostream& operator<<(std::ostream &cout, const Vector &vec) {
    
    std::cout << "{ ";
    for (const auto &el: vec) std::cout << el << ", ";
    std::cout << " }";

    return cout;
}