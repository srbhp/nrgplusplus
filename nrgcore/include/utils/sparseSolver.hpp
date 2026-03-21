#pragma once
#include "functions.hpp"
#include "qmatrix.hpp"
#include "timer.hpp"
#include <algorithm>
#include <cstddef>
#define EIGEN_USE_MKL_ALL
#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <vector>
// Spectra - header-only sparse eigenvalue solver
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>

/**
 * @brief Exact dense matrix eigenvalue solver for small systems.
 *
 * Solves the full eigenvalue problem for a dense matrix using LAPACK.
 * Suitable for small to medium systems where memory/speed trade-off favors
 * dense storage. Uses LAPACK's dgeev for non-symmetric eigenvalue problems.
 *
 * The solver finds the eigenvector corresponding to the eigenvalue with
 * minimum absolute value, normalized so that components sum to 1.
 *
 * @note Eigenvector is returned as the ground state wavefunction
 * @note Uses O(N^3) LAPACK algorithms - suitable for N < 5000
 */
class exactSolver {
  /// @brief Dense matrix storage for Hamiltonian
  qmatrix<double> qmat;
  /// @brief Matrix dimension (rows = columns, assuming square matrix)
  size_t          column;

public:
  /**
   * @brief Initialize solver with matrix size.
   *
   * @param n Dimension of square matrix (default: 0)
   */
  explicit exactSolver(size_t n = 0) : column(n) { qmat.resize(n, n, 0); }

  /**
   * @brief Add a matrix element (cumulative, supports += operation).
   *
   * @param i Row index
   * @param j Column index
   * @param val Value to add to (i,j) element
   *
   * @note This method accumulates: a(i,j) += val
   */
  void set(size_t i, size_t j, double val) {
    qmat(i, j) = qmat(i, j) + val;
  }

  /**
   * @brief Solve eigenvalue problem and return ground state eigenvector.
   *
   * Finds left/right eigenvector pairs and eigenvalues using non-symmetric
   * LAPACK routines. Returns the right eigenvector corresponding to the
   * eigenvalue with smallest absolute value, normalized to unit sum.
   *
   * @return Vector containing normalized ground state eigenvector components
   * @throw std::exception if LAPACK eigenvalue computation fails
   *
   * @note Output is normalized: sum(components) = 1.0
   * @note Modifies internal matrix during computation
   */
  std::vector<double> solve() {
    auto                diagSolver = qmat.nonsys_diag_real();
    auto                leftVec    = std::get<0>(diagSolver);
    auto                rightVec   = std::get<1>(diagSolver);
    auto                eig        = std::get<2>(diagSolver);
    std::vector<double> X(column, {0});
    double              trace = 0;
    size_t              idx   = minIndex(eig); // index of minimum absolute value
#pragma omp parallel for reduction(+ : trace)
    for (size_t i = 0; i < column; i++) {
      X[i] = rightVec(i, idx).real();
      trace += X[i];
    }
    std::cout << "EigenValue" << eig[idx] << " trace :" << trace << std::endl;
    std::transform(X.begin(), X.end(), X.begin(), [trace](double x) { return x / trace; });
    return X;
  }

  /**
   * @brief Compute matrix-vector product y = H·x.
   *
   * @param a Input vector of length column
   * @return Result vector containing H·a
   */
  auto vectorDot(std::vector<double> &a) {
    auto tmMat = qmat.dot(qmatrix<double>(a, column, 1));
    return std::vector<double>(tmMat.begin(), tmMat.end());
  }
};

/**
 * @brief Sparse matrix eigenvalue solver using Eigen + Spectra libraries.
 *
 * Solves sparse eigenvalue problems efficiently for large systems.
 * Uses Spectra library for sparse iterative eigensolvers with shift-invert
 * method to find eigenvalues near a target shift point.
 *
 * The solver finds the eigenvector with eigenvalue having minimum absolute
 * value (ground state), normalized so components sum to 1.
 *
 * @note Uses Spectra shift-invert which is efficient for interior eigenvalues
 * @note Memory-efficient for sparse matrices with O(nnz) storage
 * @note Suitable for large systems N > 5000 with low density
 * @note Requires matrix to be stored in COO (coordinate) format internally
 */
class sparseEigen {
  using datatype = double;
  /// @brief Triplet list for efficient sparse matrix assembly (COO format)
  std::vector<Eigen::Triplet<datatype>> coefficients;
  /// @brief Matrix dimension (rows = columns)
  size_t                        column;
  /// @brief Eigen sparse matrix in compressed format
  Eigen::SparseMatrix<datatype> SparseA;

public:
  /**
   * @brief Initialize sparse solver with matrix size.
   *
   * @param n Dimension of square matrix (default: 0)
   */
  explicit sparseEigen(size_t n = 0)
      : column(n), SparseA(Eigen::SparseMatrix<datatype>(n, n)) {}

  /**
   * @brief Add a matrix element (cumulative).
   *
   * Elements are stored in triplet format until getMatrix() or solve()
   * is called, at which point they're compressed into CSR format.
   *
   * @param i Row index
   * @param j Column index
   * @param value Value to add at (i,j)
   *
   * @note Accumulates: a(i,j) += value if element (i,j) added multiple times
   */
  void set(size_t i, size_t j, double value) {
    coefficients.emplace_back(i, j, value);
  }

  /**
   * @brief Get the compressed sparse matrix.
   *
   * Converts triplet format to compressed sparse row (CSR) format.
   *
   * @return Reference to compressed Eigen sparse matrix
   *
   * @note Must be called before operations if matrix was modified
   */
  auto getMatrix() {
    SparseA.setFromTriplets(coefficients.begin(), coefficients.end());
    return SparseA;
  }

  /**
   * @brief Solve sparse eigenvalue problem and return ground state eigenvector.
   *
   * Uses Spectra GenEigsRealShiftSolver with shift-invert transform.
   * Finds 1 eigenvalue with magnitude closest to zero (ground state).
   *
   * @return Normalized ground state eigenvector (sum = 1.0)
   * @throw std::runtime_error If Spectra eigenvalue computation fails
   *
   * @note Number of Lanczos vectors used: column/2
   * @note Shift point: 0 (finds eigenvalue closest to 0)
   * @note Tolerance: 1e-6 for eigenvalue convergence
   */
  auto solve() {
    SparseA.setFromTriplets(coefficients.begin(), coefficients.end());
    
    Spectra::SparseGenRealShiftSolve<double> op(SparseA);
    Spectra::GenEigsRealShiftSolver          eigs(op, 1, column / 2, 1e-6);
    
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn);
    
    std::vector<double> X(column, {0});
    if (eigs.info() == Spectra::CompInfo::Successful) {
      double           rhoTrace = 0.0;
      Eigen::VectorXcd evalues  = eigs.eigenvalues();
      auto             eigVec   = eigs.eigenvectors();
      
      for (size_t i = 0; i < column; i++) {
        X[i]     = eigVec(i, 0).real();
        rhoTrace = rhoTrace + X[i];
      }
      
      for (auto &aa : X) {
        aa = aa / rhoTrace;
      }
      std::cout << "Lowest Eigenvector : " << evalues << "Trace:  " << rhoTrace << std::endl;
    } else {
      throw std::runtime_error("Failed to evaluate Eigenvalues\n");
    }
    return X;
  }

  /**
   * @brief Compute sparse matrix-vector product y = H·x.
   *
   * Manually iterates through sparse matrix structure for matrix-vector product.
   *
   * @param X Input vector of length column
   * @return Result vector containing H·X
   *
   * @note Uses manual iteration for maximum control and efficiency
   * @note O(nnz) complexity where nnz is number of non-zeros
   */
  auto vectorDot(const std::vector<double> &X) {
    SparseA.setFromTriplets(coefficients.begin(), coefficients.end());
    std::vector<double> result(column, 0);
    
    for (int k = 0; k < SparseA.outerSize(); ++k) {
      for (Eigen::SparseMatrix<datatype>::InnerIterator it(SparseA, k); it; ++it) {
        result[it.row()] += it.value() * X[it.col()];
      }
    }
    return result;
  }
};
