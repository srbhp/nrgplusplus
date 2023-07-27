
.. _program_listing_file_nrgcore_include_nrgcore_qOperator.hpp:

Program Listing for File qOperator.hpp
======================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_nrgcore_qOperator.hpp>` (``nrgcore/include/nrgcore/qOperator.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   class qOperator {
     //
     //
     //
     //
     //
     //
     std::map<std::array<size_t, 2>, qmatrix<>> storage;
   
   public:
     qOperator() = default;
     qOperator(qmatrix<double> &opr, size_t i, size_t j) {
       std::array<size_t, 2> idx = {i, j}; // ToDo : Check for overrides
       storage[idx]              = std::move(opr);
     }
     auto getMap() { return &storage; }
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
     std::optional<qmatrix<double> *> get(size_t i, size_t j) {
       std::array<size_t, 2> idx = {i, j}; // ToDo : Check for overrides
       if (storage.find(idx) != storage.end()) {
         return &storage[idx];
       }
       // else return empty storage
       return {};
     }
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
     void clear() {
       // Clear Matrix
       for (auto &aa : storage) {
         aa.second.clear();
       }
       storage.clear();
     }
     void display(std::ostream &out) const {
       for (const auto &aa : storage) {
         out << "[" << aa.first[0] << "," << aa.first[1] << "]="
             << "|" << aa.second;
       }
     }
     friend std::ostream &operator<<(std::ostream &out, const qOperator &val) {
       out << "\n";
       val.display(out);
       return out;
     }
   };
