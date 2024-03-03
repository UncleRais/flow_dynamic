#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

#include <omp.h>

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "vtu11/vtu11.hpp"

typedef double T;
typedef std::size_t size_type;
typedef Eigen::SparseMatrix<T> matrix;
typedef Eigen::VectorXd vector;
typedef Eigen::Triplet<T> triplet;

// this struct contains all information about problem:
// 1) physical characteristics 2) computational domain sizes
// 3) mesh 4) boundary condition 5) initial distributions
 
struct problem_params {
    size_type _Nx, _Ny, _steps, _size;
    size_type _count_x, _count_y, _count_k;
    T _L, _H, _time, _hx, _hy, _tau, _nu;
    std::array<T, 4> _BC_velocity;
    std::string _format = ".vtu";
    std::string _filename;

    problem_params(size_type Nx, size_type Ny, size_type steps, T L, T H, T time, T nu, 
                   std::array<T, 4> BC_velocity, std::string filename): 
    _Nx(Nx), _Ny(Ny), _steps(steps), _L(L), _H(H), _time(time), _nu(nu), _BC_velocity{BC_velocity}, _filename(filename) {
        _hx = L / Nx; 
        _hy = H / Ny; 
        _tau = time / steps;
        _count_x = _Nx + 1;
        _count_y = _Ny + 1;
        _count_k = _count_x * _count_y;
        _size = 2 * _count_k;
    };
    ~problem_params() = default;

    //mesh indexation convertation to solution vector indexation
    size_type ij_to_k(size_type i, size_type j, size_type lower) const {
        return i + j * _count_y + lower * _count_k;
    };

    //mesh indexation convertation to matrix indexation
    std::pair<size_type, size_type> ij_to_rc(size_type template_i, size_type template_j, 
                                             size_type term_i,     size_type term_j, 
                                             size_type down_or_up, size_type psi_or_omega) const {
        // row corresponds to the 'k' of main template vertex
        size_type row = ij_to_k(template_i, template_j, down_or_up);
        // col corresponds to the 'k' of the term we want to add
        size_type col = ij_to_k(term_i, term_j, psi_or_omega);
        return std::make_pair(row, col);
    };

    void write_mesh_to_vtu(const vector& sol_prev, const vector& velocity_u, const vector& velocity_v) const {

        std::vector<T> points;
        points.reserve(3 * _size);

        //here we obtain points
        for (size_type j = 0; j < _count_y; ++j)
            for (size_type i = 0; i < _count_x; ++i) {
                points.push_back(i * _hx);
                points.push_back(j * _hy);
                points.push_back(T(0.0));
            };

        //here we obtain mesh geometry
        std::vector<vtu11::VtkIndexType> connectivity;
        std::vector<vtu11::VtkIndexType> offsets;
        std::vector<vtu11::VtkCellType> types;
        connectivity.reserve(_count_k);
        offsets.reserve(_count_k);
        types.reserve(_count_k);

        size_type offset = 0;
        size_type type = 8;
        for (size_type j = 0; j < _Ny; ++j)
            for (size_type i = 0; i < _Nx; ++i) {
                offset += 4;
                connectivity.push_back(i + j * _count_x);
                connectivity.push_back(i + 1 + j * _count_x);
                connectivity.push_back(i + (j + 1) * _count_x);
                connectivity.push_back(i + 1 + (j + 1) * _count_x);
                offsets.push_back(offset);
                types.push_back(type);
            };

        // Create small proxy mesh type
        vtu11::Vtu11UnstructuredMesh mesh { points, connectivity, offsets, types };

        std::vector<T> omega;
        std::vector<T> psi;
        std::vector<T> velocity;
        omega.resize(_count_k);
        psi.resize(_count_k);
        std::copy(sol_prev.begin(),            sol_prev.begin() + _count_k, omega.begin());
        std::copy(sol_prev.begin() + _count_k, sol_prev.end(),              psi.begin());
        for (size_type k = 0; k < _count_k; ++k) {
            velocity.push_back(velocity_u[k]);
            velocity.push_back(velocity_v[k]);
            velocity.push_back(T(0.0));
        };

        // Create tuples with (name, association, number of components) for each data set
        std::vector<vtu11::DataSetInfo> dataSetInfo
        {
            { "Omega",      vtu11::DataSetType::PointData, 1 },
            { "Psi",        vtu11::DataSetType::PointData, 1 },
            { "Velocity",   vtu11::DataSetType::PointData, 3 }
        };

        // Write data to .vtu file using Ascii format
        vtu11::writeVtu( _filename + _format, mesh, dataSetInfo, { omega, psi, velocity }, "Ascii" );
    };
};

void print_system (const matrix& A, const vector& rhs) {
    assert(A.rows() == A.cols());
    assert(A.rows() == rhs.rows());
    const size_type n = A.rows();
    std::cout << "system:" << std::endl;
    for (size_type i = 0; i < n; ++i) {
        std::cout << "{ ";
        for (size_type j = 0; j < n; ++j) {
            std::cout << std::setw(7) << std::setprecision(6) << A.coeff(i, j) << " ";
        };
        std::cout << " }  |" << std::setw(7) << std::setprecision(6) << rhs.coeff(i) << "|"<< std::endl;
    };
};

std::ostream& operator <<(std::ostream& cout, const vector& vec) {
    std::cout << "{ ";
    for (auto& el: vec) {
        std::cout << el << " ";
    };
    std::cout << " }";
    return cout;
};

void build_test_problem(std::vector<triplet>& coefficients, Eigen::VectorXd& b, int n){
    for (size_type i = 0; i < n; ++i) {
        b[i] = 0.0001;
        for (size_type j = 0; j < n; ++j) 
            if(i == j) {
                const T val = (j + 1) / 1000.0;
                coefficients.push_back(triplet(i, j, val));
            };
    };
};



