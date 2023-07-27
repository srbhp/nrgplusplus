#pragma once
#include "nrgcore/qOperator.hpp"
#include "utils/qmatrix.hpp"
#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>
class twoChannel {
  /** This class is for a  two channel (or N ) orbital with spin up and down
   * f operator.  multi-channel SIAM can made from this class.
   * This class is inherited from from `spinhalf` class. So
   * we have to implement the public  (and private functions
   * to properly set  the array's and matrices ).
   *
   *
   */
public:
  // number of f operators
  twoChannel();
  std::vector<std::vector<int>>    get_basis();
  std::vector<std::vector<double>> get_eigenvaluesQ();
  std::vector<double>              get_chi_Q();
  // protected:
  // functions
  //
  std::vector<qOperator>           f_dag_operator;
  std::vector<std::vector<double>> eigenvalues_Q;
  std::vector<double>              chi_Q;
  std::vector<std::vector<int>>    n_Q;
};
