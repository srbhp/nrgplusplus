
.. _program_listing_file_nrgcore_include_nrgcore_qsymmetry.hpp:

Program Listing for File qsymmetry.hpp
======================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_nrgcore_qsymmetry.hpp>` (``nrgcore/include/nrgcore/qsymmetry.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include <iostream>
   #include <string>
   #include <vector>
   class qsymmetry {
     // This is a base class for the symmetries of
     // the hamiltonian. This class
     size_t                   no_of_symmetry{};
     std::vector<std::string> sys_string;
   
   public:
     qsymmetry() : sys_string({"qsymmetry"}) {}
     void add_symmetry(const std::string &_id) {
       no_of_symmetry++;
       sys_string.push_back(_id);
     }
     [[nodiscard]] size_t get_symmetry_size() const { return no_of_symmetry; }
     void                 print_symmetry() const {
       std::cout << "Symmetry of this model:" << std::endl;
       for (const auto &s : sys_string) {
         std::cout << "\t " << s << std::endl;
       }
     }
   };
