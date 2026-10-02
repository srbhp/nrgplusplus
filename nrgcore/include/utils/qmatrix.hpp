#pragma once
// #include "utils/timer.hpp"
#include <cmath>
#include <complex>
#include <initializer_list>
#include <iostream>
#include <numeric>

#if defined(USE_MKL) || defined(HAVE_MKL)
#ifndef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#endif
#include <mkl.h>
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif
// #include <omp.h>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>
template <typename NType1>
std::ostream &operator<<(std::ostream &out, const std::vector<NType1> &val) {
  for (auto const &aa : val) {
    out << aa << " ";
  }
  out << "\n";
  return out;
}
template <typename NType1>
std::ostream &operator<<(std::ostream                           &out,
                         const std::vector<std::vector<NType1>> &val) {
  out << "\n";
  for (auto x : val) {
    for (auto y : x) {
      out << y << " ";
    }
    out << "\n";
  }
  return out;
}
/**
 * @brief Convenience logger used for quick debugging in numerical code.
 *
 * @tparam T Type of the value to print.
 * @param name Label shown before the value.
 * @param x Value to emit to standard output.
 */
template <typename T> void LOGGER(const std::string &name, T x) {
  std::cout << "## " << name << " : " << x << std::endl;
}
// #define logger(name) LOGGER(#name, (name))

/**
 * @class qmatrix
 * @brief A template matrix class for quantum (q) operators and general linear algebra.
 *
 * qmatrix is a row-major matrix container for NRG quantum calculations, supporting:
 * - Generic template types (double, std::complex<double>)
 * - Dense matrix operations and eigenvalue decomposition
 * - Kronecker products for tensor operations
 * - Intel MKL optimized BLAS/LAPACK routines
 *
 * The matrix is stored in flat STL vector in row-major order (LAPACK_ROW_MAJOR),
 * with element (i,j) at index i*column+j.
 *
 * @tparam T Matrix element type (default: double). Supports std::complex<double>.
 * @note All indexing is 0-based. OpenMP parallelization used for large operations.
 */
using cm_vec = std::vector<std::complex<double>>;
template <class T = double> // default is double
class qmatrix {
  /// @brief Flat vector storage of matrix elements in row-major order.
  std::vector<T> mat{0};
  /// @brief Number of rows in the matrix.
  size_t         row{0};
  /// @brief Number of columns in the matrix.
  size_t         column{0};
  /// @brief Total number of stored elements.
  size_t         dim{0};

public:
  /**
   * @brief Construct a matrix from an initializer list.
   *
   * If no explicit shape is provided, the list is interpreted as a square matrix.
   *
   * @param inputVec Row-major values to initialize the matrix with.
   * @param _row Number of rows; when zero the size is inferred as square.
   * @param _column Number of columns; when zero the size is inferred as square.
   * @throws std::runtime_error If a square matrix is implied but the flat size is not square.
   */
  qmatrix(const std::initializer_list<T> inputVec, size_t _row = 0,
          size_t _column = 0)
      : mat(inputVec) {
    if (_row == 0 && _column == 0) {
      this->row    = std::sqrt(inputVec.size());
      this->column = std::sqrt(inputVec.size());
      this->dim    = inputVec.size();
      if (this->dim != this->row * this->column) {
        throw std::runtime_error(
            "initializer_list: vector is not a square matrix! ");
      }
    } else {
      this->row    = _row;
      this->column = _column;
      this->dim    = _column * _row;
    }
    this->mat.shrink_to_fit();
  }

  /**
   * @brief Construct a matrix from an STL vector.
   *
   * @param inputVec Flat row-major values.
   * @param _row Number of rows; zero implies a square matrix.
   * @param _column Number of columns; zero implies a square matrix.
   * @throws std::runtime_error If the supplied vector is not compatible with a square matrix.
   */
  qmatrix(const std::vector<T> &inputVec, size_t _row = 0, // NOLINT
          size_t _column = 0) {
    if (_row == 0 && _column == 0) {
      this->row    = std::sqrt(inputVec.size());
      this->column = std::sqrt(inputVec.size());
      this->dim    = inputVec.size();
      if (this->dim != this->row * this->column) {
        throw std::runtime_error("vector is not a square matrix! ");
      }
    } else {
      this->row    = _row;
      this->column = _column;
      this->dim    = _column * _row;
    }
    this->mat = inputVec;
    this->mat.shrink_to_fit();
  }

  /**
   * @brief Allocate a matrix with the given shape and fill value.
   *
   * @param _row Number of rows.
   * @param _column Number of columns.
   * @param populate Fill value assigned to every entry.
   */
  qmatrix(size_t _row, size_t _column, T populate) {
    this->row    = _row;
    this->column = _column;
    this->dim    = _column * _row;
    this->mat    = std::vector<T>(_row * _column, populate);
    this->mat.shrink_to_fit();
  }

  /**
   * @brief Resize the matrix and repopulate all entries.
   *
   * @param _row Number of rows in the new matrix.
   * @param _column Number of columns in the new matrix.
   * @param populate Initial value for each new element.
   */
  void resize(size_t _row = 0, size_t _column = 0, T populate = 0) {
    this->row    = _row;
    this->column = _column;
    this->dim    = _column * _row;
    this->mat    = std::vector<T>(_row * _column, populate);
    this->mat.shrink_to_fit();
  }

  /**
   * @brief Reset the matrix to an empty state.
   */
  void clear() {
    this->row    = 0;
    this->column = 0;
    this->dim    = 0;
    this->mat.clear();
    this->mat.shrink_to_fit();
  }

  /**
   * @brief Construct a square matrix with a constant fill value.
   *
   * @param N Matrix dimension.
   * @param populate Value assigned to each element.
   */
  qmatrix(size_t N, T populate) : qmatrix(N, N, populate) {}

  /**
   * @brief Default constructor; creates an empty matrix.
   */
  qmatrix() { qmatrix(0, 0, 0); }

  /**
   * @brief Access a single flattened element by linear index.
   *
   * @param i Linear storage index.
   * @return Reference to the element at that index.
   */
  [[nodiscard]] T &operator()(size_t i) { return this->mat[i]; }

  /**
   * @brief Access a single flattened element by linear index in const context.
   *
   * @param i Linear storage index.
   * @return Copy of the value at that index.
   */
  [[nodiscard]] T operator()(size_t i) const { return this->mat[i]; }

  /**
   * @brief Access a single flattened element by linear index.
   *
   * @param i Linear storage index.
   * @return Reference to the element at that index.
   */
  [[nodiscard]] T &at(size_t i) { return this->mat[i]; }

  /**
   * @brief Access a single flattened element by linear index in const context.
   *
   * @param i Linear storage index.
   * @return Copy of the value at that index.
   */
  [[nodiscard]] T at(size_t i) const { return this->mat[i]; }

  /**
   * @brief Return the number of stored elements.
   *
   * @return The matrix capacity in flat storage, equal to row * column.
   */
  [[nodiscard]] size_t size() const { return dim; }

  /**
   * @brief Return the number of rows.
   * @return Row count.
   */
  [[nodiscard]] size_t getrow() const { return row; }

  /**
   * @brief Return the number of columns.
   * @return Column count.
   */
  [[nodiscard]] size_t getcolumn() const { return column; }

  /**
   * @brief Access an element using row and column indices.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the selected element.
   */
  [[nodiscard]] T &operator()(size_t i, size_t j) {
    return this->mat[i * column + j];
  }

  /**
   * @brief Access an element using row and column indices in const context.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Copy of the selected element.
   */
  [[nodiscard]] T operator()(size_t i, size_t j) const {
    return this->mat[i * column + j];
  }

  /**
   * @brief Access an element using row and column indices with bounds-safe semantics.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the selected element.
   */
  [[nodiscard]] T &at(size_t i, size_t j) { return this->mat[i * column + j]; }

  /**
   * @brief Access an element using row and column indices with const semantics.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Copy of the selected element.
   */
  [[nodiscard]] T at(size_t i, size_t j) const {
    return this->mat[i * column + j];
  }

  /**
   * @brief Sum all matrix elements.
   *
   * @return Sum of every entry in the matrix.
   */
  [[nodiscard]] T sum() const {
    return std::accumulate(this->mat.begin(), this->mat.end(), T{});
  }

  /**
   * @brief Sum the absolute values of all matrix entries.
   *
   * @return Absolute-value sum.
   */
  [[nodiscard]] T absSum() const {
    double sum2 = 0;
    for (auto aa : this->mat) {
      sum2 = sum2 + std::fabs(aa);
    }
    return sum2;
  }

  /**
   * @brief Compute the trace of the matrix.
   *
   * @return Sum of the diagonal elements.
   * @throws std::runtime_error If the matrix is not square.
   */
  [[nodiscard]] T trace() const {
    if (this->column == this->row) {
      T result{0};
      for (size_t i = 0; i < this->row; i++) {
        result += this->mat[i * column + i];
      }
      return result;
    }
    throw std::runtime_error("Matrix is not square matrix");
  }

  /**
   * @brief Extract the diagonal entries as a vector.
   *
   * @return Vector containing the main diagonal.
   */
  [[nodiscard]] auto getdiagonal() {
    std::vector<T> result(this->row, 0);
    for (size_t i = 0; i < this->row; i++) {
      result[i] = this->at(i, i);
    }
    return result;
  }

  /**
   * @brief Create an identity matrix of the requested size.
   *
   * @param _row Matrix dimension. When zero, the current square dimension is used.
   * @return Identity matrix with dimension _row x _row.
   * @throws std::runtime_error If the current matrix is not square.
   */
  [[nodiscard]] qmatrix<T> id(size_t _row = 0) const {
    if (this->row != this->column) {
      throw std::runtime_error("Matrix is not square matrix");
    }
    if (_row == 0) {
      _row = this->row;
    }
    qmatrix<T> result(_row, _row, 0);
    for (size_t i = 0; i < _row; i++) {
      result(i, i) = 1.0;
    }
    return result;
  }

  /**
   * @brief Extract the real part of a complex-valued matrix.
   *
   * @return Matrix containing real parts of each entry.
   */
  [[nodiscard]] qmatrix<double> real() const {
    qmatrix<double> result(this->column, this->row, 0);
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->at(i).real();
    }
    return result;
  }

  /**
   * @brief Extract the imaginary part of a complex-valued matrix.
   *
   * @return Matrix containing imaginary parts of each entry.
   */
  [[nodiscard]] qmatrix<double> imag() const {
    qmatrix<double> result(this->column, this->row, 0);
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->at(i).imag();
    }
    return result;
  }

  /**
   * @brief Compute the conjugate transpose of the matrix.
   *
   * @return Transposed matrix with complex conjugation applied.
   */
  [[nodiscard]] qmatrix<T> cTranspose() const {
    qmatrix<T> result(this->column, this->row, 0);
    if constexpr (std::is_same_v<T, std::complex<double>>) {
      for (size_t i = 0; i < this->row; i++) {
        for (size_t j = 0; j < this->column; j++) {
          result(j, i) = std::conj(this->mat[i * column + j]);
        }
      }
    } else {
      for (size_t i = 0; i < this->row; i++) {
        for (size_t j = 0; j < this->column; j++) {
          result(j, i) = this->mat[i * column + j];
        }
      }
    }
    return result;
  }

  /**
   * @brief Multiply the matrix by a scalar.
   *
   * @param x Scalar multiplier.
   * @return Result of element-wise scalar multiplication.
   */
  [[nodiscard]] qmatrix operator*(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) * x;
    }
    return result;
  }

  /**
   * @brief Divide the matrix by a scalar.
   *
   * @param x Scalar divisor.
   * @return Result of element-wise scalar division.
   */
  [[nodiscard]] qmatrix operator/(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) / x;
    }
    return result;
  }

  /**
   * @brief Return a pointer to the underlying flat storage.
   *
   * @return Pointer to the internal data array.
   */
  T *data() { return this->mat.data(); }

  /**
   * @brief Return a const pointer to the underlying flat storage.
   *
   * @return Const pointer to the internal data array.
   */
  [[nodiscard]] const T *data() const { return this->mat.data(); }

  /**
   * @brief Return a const iterator to the beginning of the underlying storage.
   * @return Begin iterator.
   */
  [[nodiscard]] auto begin() const { return this->mat.begin(); }

  /**
   * @brief Return a const iterator to the end of the underlying storage.
   * @return End iterator.
   */
  [[nodiscard]] auto end() const { return this->mat.end(); }

  /**
   * @brief Print the matrix contents to standard output.
   */
  void display() {}

  /**
   * @brief Add two matrices of identical shape.
   *
   * @param rhs Matrix on the right-hand side.
   * @return Element-wise sum.
   * @throws std::runtime_error If the matrix dimensions differ.
   */
  [[nodiscard]] qmatrix operator+(const qmatrix<T> &rhs) const {
    if (this->row == rhs.row && this->column == rhs.column) {
      qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
      for (size_t i = 0; i < this->dim; i++) {
        result(i) = this->at(i) + rhs(i);
      }
      return result;
    }
    throw std::runtime_error("qmatrix have different size for operator+");
  }

  /**
   * @brief Subtract two matrices of identical shape.
   *
   * @param rhs Matrix on the right-hand side.
   * @return Element-wise difference.
   * @throws std::runtime_error If the matrix dimensions differ.
   */
  [[nodiscard]] qmatrix operator-(const qmatrix<T> &rhs) const {
    if (this->row == rhs.row && this->column == rhs.column) {
      qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
      for (size_t i = 0; i < this->dim; i++) {
        result(i) = this->mat[i] - rhs(i);
      }
      return result;
    }
    throw std::runtime_error("qmatrix have different size for operator -\n");
  }

  /**
   * @brief Multiply this matrix by another matrix using BLAS-backed matrix multiplication.
   *
   * @param rhs Right-hand matrix with compatible dimensions.
   * @param talpha Optional scalar prefactor.
   * @return Product matrix of dimension row x rhs.column.
   * @throws std::runtime_error If the inner dimensions do not match.
   */
  [[nodiscard]] qmatrix<T> dot(const qmatrix<T> &rhs, double talpha = 1.0) {
    if (this->column == rhs.row) {
      qmatrix result(this->row, rhs.column, 0);
      size_t  m = this->row;
      size_t  k = this->column;
      size_t  n = rhs.column;
      if (m == 0 || k == 0 || n == 0) {
        return result;
      }
      if constexpr (std::is_same_v<T, double>) {
        const double alpha = talpha;
        const double beta  = 0;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                    this->data(), k, rhs.data(), n, beta, result.data(), n);
      }
      if constexpr (std::is_same_v<T, std::complex<double>>) {
        const std::complex<double> alpha{talpha, 0.0};
        const std::complex<double> beta{0.0, 0.0};
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
                    this->data(), k, rhs.data(), n, &beta, result.data(), n);
      }
      return result;
    }
    throw std::runtime_error("dot:qmatrix have different size for dot\n");
  }

  /**
   * @brief Diagonalize a symmetric or Hermitian matrix and return eigenvalues.
   *
   * @return Vector of eigenvalues in ascending order as returned by LAPACK.
   * @throws std::runtime_error If the matrix is not square.
   */
  [[nodiscard]] std::vector<double> diag() {
    if (this->row != this->column) {
      throw std::runtime_error("Error: Matrix is not a square matrix! \n");
    }
    std::vector<double> w(this->row, 0);
    size_t              n = w.size();
    if (n == 0) {
      return w;
    }
    int info = -1;
    if constexpr (std::is_same_v<T, double>) {
      info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, this->mat.data(), n,
                            w.data());
    }
    if constexpr (std::is_same_v<T, std::complex<double>>) {
      info = LAPACKE_zheevd(
          LAPACK_ROW_MAJOR, 'V', 'U', n,
          this->mat.data(), n, w.data());
    }
    if (info > 0) {
      std::cout << "Error:Not able to solve Eigen value problem." << std::endl;
    }
    return w;
  }

  /**
   * @brief Solve the nonsymmetric eigenproblem for a complex matrix.
   *
   * @return Tuple containing left eigenvectors, right eigenvectors, and eigenvalues.
   * @throws std::runtime_error If the matrix is not square.
   * @throws std::invalid_argument If the matrix is not complex-valued.
   */
  [[nodiscard]] std::tuple<qmatrix<T>, qmatrix<T>, cm_vec>
  nonsys_diag_complex() {
    if (this->row != this->column) {
      throw std::runtime_error("Error: Matrix is not a square matrix! \n");
    }
    std::vector<T> w(this->row, 0);
    size_t         n = w.size();
    qmatrix<T>     lv(n, n, 0);
    qmatrix<T>     rv(n, n, 0);
    if constexpr (std::is_same_v<T, std::complex<double>>) {
      auto info = LAPACKE_zgeev(
          LAPACK_ROW_MAJOR, 'V', 'V', n,
          this->data(), n,
          w.data(),
          lv.data(), n,
          rv.data(), n);
#pragma omp parallel for // NOLINT
      for (size_t i = 0; i < n; i++) {
        std::complex<double> aa{0};
        for (size_t k = 0; k < n; k++) {
          aa += std::conj(lv(k, i)) * rv(k, i);
        }
        for (size_t k = 0; k < n; k++) {
          lv(k, i) = lv(k, i) / std::conj(std::sqrt(aa));
          rv(k, i) = rv(k, i) / (std::sqrt(aa));
        }
      }
      if (info > 0) {
        throw std::runtime_error(
            "The algorithm LAPACKE_zgeev failed to compute eigenvalues.\n");
      }
    } else {
      throw std::invalid_argument(
          "nonsys_diag_complex: This function is for complex matrices");
    }
    return std::tuple(lv, rv, w);
  }

  /**
   * @brief Solve the nonsymmetric real eigenproblem for a real-valued matrix.
   *
   * @return Tuple containing left vectors, right vectors, and complex eigenvalues.
   * @throws std::invalid_argument If the matrix is not square or is not a qmatrix<double>.
   */
  std::tuple<qmatrix<std::complex<T>>, qmatrix<std::complex<T>>,
             std::vector<std::complex<T>>>
  nonsys_diag_real() {
    size_t nsize = this->getrow();
    if (this->getrow() != this->getcolumn()) {
      throw std::invalid_argument(
          "nonsys_diag_real: This is not a square matrix");
    }
    if constexpr (!std::is_same_v<T, double>) {
      throw std::invalid_argument(
          "nonsys_diag_real: This is  not a qmatrix<double>! ");
    }
    std::vector<T> wr(nsize, 0);
    std::vector<T> wi(nsize, 0);
    std::vector<T> vl(nsize * nsize, 0);
    std::vector<T> vr(nsize * nsize, 0);
    auto info =
        LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', nsize, this->data(), nsize,
                      wr.data(), wi.data(), vl.data(), nsize, vr.data(), nsize);
    if (info > 0) {
      std::cout << "The algorithm failed to compute eigenvalues." << std::endl;
      exit(1);
    }
    std::vector<std::complex<T>> eigenvalues(nsize, 0);
    qmatrix<std::complex<T>>     leftVectors(nsize, nsize, 0);
    qmatrix<std::complex<T>>     rightVectors(nsize, nsize, 0);
#pragma omp parallel for // NOLINT
    for (size_t j = 0; j < nsize; j++) {
      eigenvalues[j] = std::complex<T>(wr[j], wi[j]);
    }
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < nsize; i++) {
      size_t j = 0;
      while (j < nsize) {
        if (wi[j] == static_cast<T>(0.0)) {
          leftVectors(i, j)  = vl[i * nsize + j];
          rightVectors(i, j) = vr[i * nsize + j];
          j++;
        } else {
          leftVectors(i, j) =
              std::complex<T>(vl[i * nsize + j], vl[i * nsize + j + 1]);
          leftVectors(i, j + 1) =
              std::complex<T>(vl[i * nsize + j], -vl[i * nsize + j + 1]);
          rightVectors(i, j) =
              std::complex<T>(vr[i * nsize + j], vr[i * nsize + j + 1]);
          rightVectors(i, j + 1) =
              std::complex<T>(vr[i * nsize + j], -vr[i * nsize + j + 1]);
          j += 2;
        }
      }
    }
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < nsize; i++) {
      std::complex<T> aa{0};
      for (size_t k = 0; k < nsize; k++) {
        aa += std::conj(leftVectors(k, i)) * rightVectors(k, i);
      }
      if (std::fabs(aa) < 1e-5) {
      } else {
        for (size_t k = 0; k < nsize; k++) {
          leftVectors(k, i)  = leftVectors(k, i) / std::conj(std::sqrt(aa));
          rightVectors(k, i) = rightVectors(k, i) / (std::sqrt(aa));
        }
      }
    }
    return {leftVectors, rightVectors, eigenvalues};
  }

  /**
   * @brief Stream a matrix to an output stream for debugging and logging.
   *
   * @param out Output stream.
   * @param val Matrix to print.
   * @return Reference to the same output stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const qmatrix<T> &val) {
    out << "\n";
    for (size_t i = 0; i < val.row; ++i) {
      for (size_t j = 0; j < val.column; ++j) {
        out << val(i, j) << " ";
      }
      out << "\n";
    }
    return out;
  }

  /**
   * @brief Compute the Kronecker product between this matrix and another.
   *
   * @param rhs Right-hand matrix.
   * @param alpha Optional scalar prefactor.
   * @return Kronecker product matrix of dimension (row * rhs.row) by (column * rhs.column).
   */
  qmatrix<T> krDot(const qmatrix<T> &rhs, double alpha = 1) {
    size_t     m = this->row;
    size_t     n = this->column;
    size_t     p = rhs.row;
    size_t     q = rhs.column;
    qmatrix<T> result(this->row * rhs.row, this->column * rhs.column, 0);
#pragma omp parallel for // NOLINT
    for (size_t r = 0; r < m; r++) {
      for (size_t s = 0; s < n; s++) {
        for (size_t v = 0; v < p; v++) {
          for (size_t w = 0; w < q; w++) {
            result((r * p + v), (s * q + w)) =
                this->at(r, s) * rhs(v, w) * alpha;
          }
        }
      }
    }
    return result;
  }

  /**
   * @brief Apply a unitary similarity transform using a supplied eigenvector basis.
   *
   * @param eigen_vector Matrix containing the eigenvectors to apply.
   */
  void unitary_transform(const qmatrix<T> &eigen_vector) {
    auto result = eigen_vector.cTranspose().dot(this->dot(eigen_vector));
    *this       = result;
  }
};
