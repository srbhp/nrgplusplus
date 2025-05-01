#pragma once
// #include "utils/timer.hpp"
#include <cmath>
#include <complex>
#include <initializer_list>
#include <iostream>
#include <numeric>
// #include <lapacke.h> // comment this if mkl is used
// #include <cblas.h>
// uncomment if you want to use the mkl class
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#include <mkl_lapacke.h>
// typedef double _Complex DCOMPLEX
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
template <typename T> void LOGGER(const std::string &name, T x) {
  std::cout << "## " << name << " : " << x << std::endl;
}
// #define logger(name) LOGGER(#name, (name))
// This is a matrix class quamtum(q) operators
// This is also consistent with Lapack_row_MAJOR
using cm_vec = std::vector<std::complex<double>>;
template <class T = double> // default is double
class qmatrix {
  std::vector<T> mat{0};
  size_t         row{0}, column{0}, dim{0};

public:
  qmatrix(const std::initializer_list<T> inputVec, size_t _row = 0,
          size_t _column = 0)
      : mat(inputVec) {
    // check if this is a quare matrix
    // std::cout << "constructed with a " << inputVec.size() << "-element
    // list\n";
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
  qmatrix(const std::vector<T> &inputVec, size_t _row = 0, // NOLINT
          size_t _column = 0) {
    // check if this is a quare matrix
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
  qmatrix(size_t _row, size_t _column, T populate) {
    this->row    = _row;
    this->column = _column;
    this->dim    = _column * _row;
    this->mat    = std::vector<T>(_row * _column, populate);
    this->mat.shrink_to_fit();
    // TODO(sp): restrict memory usage
  }
  void resize(size_t _row = 0, size_t _column = 0, T populate = 0) {
    this->row    = _row;
    this->column = _column;
    this->dim    = _column * _row;
    this->mat    = std::vector<T>(_row * _column, populate);
    this->mat.shrink_to_fit();
    // TODO(sp): restrict memory usage
  }
  void clear() {
    this->row    = 0;
    this->column = 0;
    this->dim    = 0;
    this->mat.clear();
    this->mat.shrink_to_fit();
  } // This is for square matrix
  qmatrix(size_t N, T populate) { qmatrix(N, N, populate); }
  qmatrix() { qmatrix(0, 0, 0); }
  [[nodiscard]] T     &operator()(size_t i) { return this->mat[i]; }
  [[nodiscard]] T      operator()(size_t i) const { return this->mat[i]; }
  [[nodiscard]] T     &at(size_t i) { return this->mat[i]; }
  [[nodiscard]] T      at(size_t i) const { return this->mat[i]; }
  [[nodiscard]] size_t size() const { return dim; }
  [[nodiscard]] size_t getrow() const { return row; }
  [[nodiscard]] size_t getcolumn() const { return column; }
  [[nodiscard]] T     &operator()(size_t i, size_t j) {
        return this->mat[i * column + j];
  }
  [[nodiscard]] T operator()(size_t i, size_t j) const {
    return this->mat[i * column + j];
  }
  // add a at operator
  [[nodiscard]] T &at(size_t i, size_t j) { return this->mat[i * column + j]; }
  [[nodiscard]] T  at(size_t i, size_t j) const {
     return this->mat[i * column + j];
  }
  // TODO(sp): arithmetic operator
  // TODO(sp): use stl
  [[nodiscard]] T sum() const {
    // Sum of all elements
    return std::accumulate(this->mat.begin(), this->mat.end(), 0);
  }
  [[nodiscard]] T absSum() const {
    // Sum of all absolute values of the elements
    double sum2 = 0;
    for (auto aa : this->mat) {
      sum2 = sum2 + std::fabs(aa);
    }
    return sum;
  }
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
  [[nodiscard]] auto getdiagonal() {
    std::vector<T> result(this->row, 0);
    for (size_t i = 0; i < this->row; i++) {
      result[i] = this->at(i, i);
    }
    return result;
  }
  [[nodiscard]] qmatrix<T> id(size_t _row = 0) const {
    if (this->row != this->column) {
      throw std::runtime_error("Matrix is not square matrix");
    }
    if (_row == 0) {
      _row = this->row;
    }
    qmatrix<T> result(_row, _row, 0); // reverse the row and column
    for (size_t i = 0; i < _row; i++) {
      result(i, i) = 1.0;
    }
    return result;
  }
  [[nodiscard]] qmatrix<double> real() const {
    qmatrix<double> result(this->column, this->row,
                           0); // reverse the row and column
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->at(i).real();
    }
    return result;
  }
  [[nodiscard]] qmatrix<double> imag() const {
    qmatrix<double> result(this->column, this->row,
                           0); // reverse the row and column
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->at(i).imag();
    }
    return result;
  }
  [[nodiscard]] qmatrix<T> cTranspose() const {    // Conjugate transpose
    qmatrix<T> result(this->column, this->row, 0); // reverse the row and column
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
  [[nodiscard]] qmatrix operator*(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) * x;
    }
    return result;
  }
  [[nodiscard]] qmatrix operator/(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) / x;
    }
    return result;
  }
  //  void reset(const T &x = 0) {
  //    for (auto &aa : this->mat) {
  //      aa = x;
  //    }
  //  }
  T                     *data() { return this->mat.data(); }
  [[nodiscard]] const T *data() const { return this->mat.data(); }
  [[nodiscard]] auto     begin() const { return this->mat.begin(); }
  [[nodiscard]] auto     end() const { return this->mat.end(); }
  //
  //
  void                  display() {}
  [[nodiscard]] qmatrix operator+(const qmatrix<T> &rhs) const {
    // std::cout << "Started operator+";
    if (this->row == rhs.row && this->column == rhs.column) {
      qmatrix result(this->row, this->column, 0);
#pragma omp parallel for // NOLINT
      for (size_t i = 0; i < this->dim; i++) {
        result(i) = this->at(i) + rhs(i);
      }
      // std::cout << "End of operator+" << std::endl;
      return result;
    } //
      // else throw
    throw std::runtime_error("qmatrix have different size for operator+");
  }
  // Subtract
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
  [[nodiscard]] qmatrix<T> dot(const qmatrix<T> &rhs, double talpha = 1.0) {
    if (this->column == rhs.row) {
      qmatrix result(this->row, rhs.column, 0);
      size_t  m = this->row;
      size_t  k = this->column;
      size_t  n = rhs.column;
      if (m == 0 || k == 0 || n == 0) {
        // std::cout << "One of the dimension is zero " << std::endl;
        return result;
      }
      const double alpha = talpha;
      const double beta  = 0;
      if constexpr (std::is_same_v<T, double>) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                    this->data(), k, rhs.data(), n, beta, result.data(), n);
      }
      if constexpr (std::is_same_v<T, std::complex<double>>) {
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
                    this->data(), k, rhs.data(), n, &beta, result.data(), n);
      }
      return result;
    }
    throw std::runtime_error("dot:qmatrix have different size for dot\n");
  }
  [[nodiscard]] std::vector<double> diag() {
    // This function diagonalizes a symmetric/harmitian matrix.
    // So eigen value are always real(double).
    // If the matrix is not symmetric the call nonsys_diag
    if (this->row != this->column) {
      throw std::runtime_error("Error: Matrix is not a square matrix! \n");
    }
    std::vector<double> w(this->row, 0);
    size_t              n = w.size();
    if (n == 0) {
      // std::cout << "This is a empty matrix" << std::endl;
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
          // reinterpret_cast<__complex__ double *>(this->mat.data()), n,
          this->mat.data(), n, w.data());
    }
    // int info= LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, n, w );
    if (info > 0) {
      std::cout << "Error:Not able to solve Eigen value problem." << std::endl;
    }
    return w;
  }
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
          // recast the complex pointer
          // reinterpret_cast<__complex__ double *>(this->data()), n,
          this->data(), n,
          // recast the complex pointer
          w.data(),
          // recast the complex pointer
          lv.data(), n,
          // recast the complex pointer
          rv.data(), n);
/* Check for convergence */
// Normalize the vectors only for the diagonal elements
#pragma omp parallel for // NOLINT
      for (size_t i = 0; i < n; i++) {
        std::complex<double> aa{0};
        for (size_t k = 0; k < n; k++) {
          aa += std::conj(lv(k, i)) * rv(k, i);
        }
        // if (std::fabs(aa) < 1e-5) {
        //  // std::cout << "Warning:Normed:" << aa;
        //} else {
        for (size_t k = 0; k < n; k++) {
          lv(k, i) = lv(k, i) / std::conj(std::sqrt(aa));
          rv(k, i) = rv(k, i) / (std::sqrt(aa));
        }
        //}
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
    // return {lv, rv, w};
  }
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
    //
    // timer t1("Solving exact");
    // mkl_set_num_threads(200);
    auto info =
        LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', nsize, this->data(), nsize,
                      wr.data(), wi.data(), vl.data(), nsize, vr.data(), nsize);
    // std::cout << "ExactSolver Done " << t1.getDuration() << std::endl;
    if (info > 0) {
      std::cout << "The algorithm failed to compute eigenvalues." << std::endl;
      exit(1);
    }
    std::vector<std::complex<T>> eigenvalues(nsize, 0);
    qmatrix<std::complex<T>>     leftVectors(nsize, nsize, 0);
    qmatrix<std::complex<T>>     rightVectors(nsize, nsize, 0);
// set the values
#pragma omp parallel for // NOLINT
    for (size_t j = 0; j < nsize; j++) {
      eigenvalues[j] = std::complex<T>(wr[j], wi[j]);
    }
// set eigenvector
// I dont know why the fuck this is organized this way
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
// NOrmalize the vector
#pragma omp parallel for // NOLINT
    for (size_t i = 0; i < nsize; i++) {
      std::complex<T> aa{0};
      for (size_t k = 0; k < nsize; k++) {
        aa += std::conj(leftVectors(k, i)) * rightVectors(k, i);
      }
      if (std::fabs(aa) < 1e-5) {
        // std::cout << "Warning:Normed:" << std::fabs(aa);
      } else {
        for (size_t k = 0; k < nsize; k++) {
          leftVectors(k, i)  = leftVectors(k, i) / std::conj(std::sqrt(aa));
          rightVectors(k, i) = rightVectors(k, i) / (std::sqrt(aa));
          //  aa2 += std::conj(lv(k, i)) * rv(k, j);
        }
      }
      //       std::cout << std::endl;
    }
    return {leftVectors, rightVectors, eigenvalues};
  }
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
  qmatrix<T> krDot(const qmatrix<T> &rhs, double alpha = 1) {
    // https://en.wikipedia.org/wiki/Kronecker_product
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
  void unitary_transform(const qmatrix<T> &eigen_vector) {
    // U^T. x . U
    // TODO(sp):
    auto result = eigen_vector.cTranspose().dot(this->dot(eigen_vector));
    *this       = result;
  }
};
