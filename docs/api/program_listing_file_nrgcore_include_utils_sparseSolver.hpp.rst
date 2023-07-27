
.. _program_listing_file_nrgcore_include_utils_sparseSolver.hpp:

Program Listing for File sparseSolver.hpp
=========================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_utils_sparseSolver.hpp>` (``nrgcore/include/utils/sparseSolver.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   // Spcetra
   #include <Spectra/GenEigsRealShiftSolver.h>
   #include <Spectra/MatOp/SparseGenRealShiftSolve.h>
   class exactSolver {
     qmatrix<double> qmat;
     size_t          column;
   
   public:
     explicit exactSolver(size_t n = 0) : column(n) { qmat.resize(n, n, 0); }
     void set(size_t i, size_t j, double val) {
       qmat(i, j) = qmat(i, j) + val; // +=
     }
     std::vector<double> solve() {
       // std::cout << "Solving exact" << std::endl;
       // auto [leftVec, rightVec, eig] = qmat.nonsys_diag_real();
       // Clang complains about this line, but it works.
       // auto [leftVec, rightVec, eig] = qmat.nonsys_diag_real();
       auto                diagSolver = qmat.nonsys_diag_real();
       auto                leftVec    = std::get<0>(diagSolver);
       auto                rightVec   = std::get<1>(diagSolver);
       auto                eig        = std::get<2>(diagSolver);
       std::vector<double> X(column, {0});
       double              trace = 0;
       size_t              idx   = minIndex(eig); // index of minimum value
   #pragma omp parallel for reduction(+ : trace)
       for (size_t i = 0; i < column; i++) {
         // std::cout << "Left Eigen vaector " << leftVec(i, idx) << std::endl;
         X[i] = rightVec(i, idx).real();
         trace += X[i];
       }
       std::cout << "EigenValue" << eig[idx] << " trace :" << trace << std::endl;
       // std::cout << "Evec:" << X << std::endl;
       // trace should be one
       // parallelization
       std::transform(
           // std::execution::par_unseq,
           X.begin(), X.end(), X.begin(), [trace](double x) { return x / trace; });
       return X;
     }
     auto vectorDot(std::vector<double> &a) {
       // std::cout << "Started vectorDot" << std::endl;
       auto tmMat = qmat.dot(qmatrix<double>(a, column, 1));
       // std::cout << "End of vectorDot" << std::endl;
       return std::vector<double>(tmMat.begin(), tmMat.end());
     }
   };
   /*
   class sparseArma {
     arma::sp_mat sparseMat;
     size_t       column; // Assumes that the row = column
   public:
     sparseArma(size_t n = 0) {
       column    = n;
       sparseMat = arma::sp_mat(n, n);
     }
     void set(size_t i, size_t j, double val) { sparseMat(j, i) = val; }
     auto getMatrix() { return sparseMat; }
     auto solve() {
       arma::cx_vec    eigval;
       arma::cx_mat    eigvec;
       arma::eigs_opts opts;
       opts.maxiter = 10000; // increase max iterations to 10000
       eigs_gen(eigval, eigvec, sparseMat, 1, 1e-8, opts);
       std::vector<double> X(column, {0});
       double              trace = 0;
       for (size_t i = 0; i < column; i++) {
         X[i] = eigvec(i, 0).real();
         trace += X[i];
       }
       // trace should be one
       for (auto &aa : X) {
         aa = aa / trace;
       }
       std::cout << "Eigenvalues: " << eigval << " trace " << trace << std::endl;
       return X;
     }
     auto vectorDot(const std::vector<double> &X) {
       std::cout << "Started vectorDot" << std::endl;
       auto result = arma::vec(sparseMat * arma::vec(X));
       std::cout << "Type: " << typeid(result).name() << std::endl;
       std::cout << "Endend vectorDot" << std::endl;
       return result;
     }
   };
   
   */
   class sparseEigen {
     using datatype = double;
     std::vector<Eigen::Triplet<datatype>>
                                   coefficients; // coefficients of the matrix
     size_t                        column;
     Eigen::SparseMatrix<datatype> SparseA;
   
   public:
     explicit sparseEigen(size_t n = 0)
         : column(n), SparseA(Eigen::SparseMatrix<datatype>(n, n)) {
       // column = n; // Assumes that the row = column
     }
     void set(size_t i, size_t j, double value) {
       coefficients.emplace_back(i, j, value);
     }
     auto getMatrix() {
       SparseA.setFromTriplets(coefficients.begin(), coefficients.end());
       return SparseA;
     }
     auto solve() {
       // set the matrix
       SparseA.setFromTriplets(coefficients.begin(), coefficients.end());
       // Setup solver
       // Spectra::SparseGenMatProd<datatype> op(SparseA);
       // Spectra::GenEigsSolver              eigs(op, 1, column / 2);
       Spectra::SparseGenRealShiftSolve<double> op(SparseA);
       Spectra::GenEigsRealShiftSolver          eigs(op, 1, column / 2,
                                                     1e-6); // Initialize and compute
       eigs.init();
       eigs.compute(
           Spectra::SortRule::LargestMagn); // returns the smallest eigenvalues
       // std::cout << "ncov : " << nconv << std::endl;
       // Retrieve results
       std::vector<double> X(column, {0});
       if (eigs.info() == Spectra::CompInfo::Successful) {
         double           rhoTrace = 0.0;
         Eigen::VectorXcd evalues  = eigs.eigenvalues();
         auto             eigVec   = eigs.eigenvectors(); // This change this shit
         // std::cout<< "evalues : " <<  eigVec.size() <<" "<< column<< std::endl;
         for (size_t i = 0; i < column; i++) {
           X[i]     = eigVec(i, 0).real();
           rhoTrace = rhoTrace + X[i];
         }
         // trace should be one
         for (auto &aa : X) {
           aa = aa / rhoTrace;
         }
         std::cout << "Lowest Eigenvector : " << evalues << "Trace:  " << rhoTrace
                   << std::endl;
       } else {
         throw std::runtime_error("Failed to evaluate Eigenvaleus\n");
       }
       return X;
       // return eigs.eigenvectors()
       //
     }
     auto vectorDot(const std::vector<double> &X) {
       // timer t1("vectorDot Sparse Solver");
       SparseA.setFromTriplets(coefficients.begin(), coefficients.end());
       std::vector<double> result(column, 0);
       for (int k = 0; k < SparseA.outerSize(); ++k) {
         for (Eigen::SparseMatrix<datatype>::InnerIterator it(SparseA, k); it;
              ++it) {
           // std::cout << "(" << it.row() << ","; // row index
           // std::cout << it.col() << ")\t"; // col index (here it is equal to k)
           result[it.row()] += it.value() * X[it.col()];
         }
       }
       return result;
     }
   };
