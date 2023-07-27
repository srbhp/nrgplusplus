
.. _program_listing_file_nrgcore_include_models_twoChannel.hpp:

Program Listing for File twoChannel.hpp
=======================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_models_twoChannel.hpp>` (``nrgcore/include/models/twoChannel.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
