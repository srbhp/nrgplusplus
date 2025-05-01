#pragma once
#include "utils/qmatrix.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <optional>
#include <ostream>
#include <set>
/**
 * @class qOperator
 * @brief Represents a quantum operator in a block-diagonal form.
 *
 * This class provides methods to manipulate and display quantum operators
 * used in the NRG calculations.
 */
class qOperator {
  std::map<std::array<size_t, 2>, qmatrix<>> storage;

public:
  qOperator() = default;
  qOperator(qmatrix<double> &opr, size_t i, size_t j) {
    std::array<size_t, 2> idx = {i, j}; // ToDo : Check for overrides
    storage[idx]              = std::move(opr);
  }
  /**
   * @brief Get the internal Storage object of the `qmatrix`
   *
   * @return std::map<std::array<size_t, 2>, qmatrix<>>
   */
  auto getMap() { return &storage; }
  /**
   * @brief Set the `qOperator` for the `i`th and `j`th
   * symmetry basis.
   *
   * @param opr : Matrix object
   * @param i : `i`th symmetry block
   * @param j : `j`th symmetry block
   */
  void set(const qmatrix<double> &opr, size_t i, size_t j) {
    std::array<size_t, 2> idx = {i, j}; // ToDo : Check for overrides
    if (storage.find(idx) != storage.end()) {
      throw std::runtime_error("Found operator already");
    }
    storage[idx] = opr;
  }
  void set(qmatrix<double> &opr, size_t i, size_t j) {
    // #pragma omp critical
    {
      std::array<size_t, 2> idx = {i, j}; // ToDo : Check for overrides
      if (storage.find(idx) != storage.end()) {
        throw std::runtime_error("Found operator already");
      }
      storage[idx] = std::move(opr);
    }
  }
  /**
   * @brief This function returns the `qmatrix` pointer for the sub-block
   * `i` and `j`. If the sub-block is not found, it returns an empty
   * `std::optional`.
   *
   * @param i The index of the block to be returned.
   * @param j The index of the block to be returned.
   * @return std::optional<qmatrix>
   */
  std::optional<qmatrix<double> *> get(size_t i, size_t j) {
    std::array<size_t, 2> idx = {i, j}; // ToDo : Check for overrides
    if (storage.find(idx) != storage.end()) {
      return &storage[idx];
    }
    // else return empty storage
    return {};
  }
  /**
   * @brief  This function unitary transform the `qOperator` by the `U` matrix.
   * The `U` matrix is a block-diagonal matrix of `qmatrix` type.
   * The `qOperator` is transformed as \f U^\dagger \cdot  qOperator \cdot U \f
   * .
   *
   * @param U The unitary matrix to transform the `qOperator`. Usually
   * we get this matrix after diagonalization of the Hamiltonian.
   */
  void unitaryTransform(qOperator &U) {
    // TODO(sp): Check if this is correct
    for (auto &aa : storage) {
      // out << "[" << aa.first[0] << "," << aa.first[1] << "]="
      //     << "|" << aa.second;
      //
      auto i    = aa.first[0];
      auto j    = aa.first[1];
      aa.second = U.get(i, i).value()->cTranspose().dot(aa.second).dot(
          *(U.get(j, j)).value());
    }
  }
  // std::optional<qmatrix<>> get_transposed_operator(size_t i, size_t j) {
  //  // TODO: We are returning the full storage
  //  auto idx = i + nBlocks * j;
  //  if (storage.find(idx) != storage.end()) {
  //    return storage.at(idx).transpose();
  //  } else {
  //    return {};
  //  }
  //}
  /**
   * @brief Clears the storage for the operator.
   *
   * This method removes all stored matrix elements and resets the operator.
   */
  void clear() {
    // Clear Matrix
    for (auto &aa : storage) {
      aa.second.clear();
    }
    storage.clear();
  }
  /**
   * @brief Prints the `qOperator` to the output stream.
   *
   * Avoid calling this function if the size of the `qOperator` is large.
   *
   * @param out Output stream to print to.
   */
  void display(std::ostream &out) const {
    for (const auto &aa : storage) {
      out << "[" << aa.first[0] << "," << aa.first[1] << "]="
          << "|" << aa.second;
    }
  }
  /**
   * @brief Overloads the stream insertion operator for `qOperator`.
   *
   * Allows the use of `std::cout << qOperator` or any other `ostream` object.
   *
   * @param out The output stream.
   * @param val The `qOperator` to print.
   * @return The modified output stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const qOperator &val) {
    out << "\n";
    val.display(out);
    return out;
  }
};
