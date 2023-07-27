
.. _program_listing_file_nrgcore_include_nrgcore_nrgData.hpp:

Program Listing for File nrgData.hpp
====================================

|exhale_lsh| :ref:`Return to documentation for file <file_nrgcore_include_nrgcore_nrgData.hpp>` (``nrgcore/include/nrgcore/nrgData.hpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "nrgcore/qOperator.hpp"
   #include "utils/h5stream.hpp"
   #include "utils/qmatrix.hpp"
   #include "utils/timer.hpp"
   #include <array>
   #include <ctime>
   #include <iostream>
   #include <map>
   #include <numeric>
   #include <optional>
   #include <random>
   #include <string>
   #include <tuple>
   #include <vector>
   template <typename nrgcore_type> class NrgData {
     std::string        tmpNrgFilename;
     nrgcore_type      *nrg_object;
     h5stream::h5stream pfiletmp;
     bool               isClosed = false;
   
   public:
     explicit NrgData(const std::string &tfilename =
                          "") { // We don't want to take nrgObject here
       if (tfilename.empty()) {
         tmpNrgFilename = "tempfile" + std::to_string(std::rand()) + ".h5";
       } else {
         tmpNrgFilename = tfilename;
       }
       pfiletmp.setFileName(tmpNrgFilename);
       // TODO(sp): clear the input operator
       nrg_object = nullptr;
     }
     void readFromFile(const std::string &tfile) {
       remove(tmpNrgFilename.c_str());
       tmpNrgFilename = tfile;
       pfiletmp.close();
       pfiletmp.setFileName(tmpNrgFilename, "r");
       pfiletmp.read<size_t>(savedNRGIndex, "savedNRGIndex");
       if (nrg_object != nullptr) {
         pfiletmp.read(nrg_object->relativeGroundStateEnergy,
                       "relativeGroundStateEnergy");
       }
     }
     void setFileName(const std::string &tfile) {
       tmpNrgFilename = tfile;
       pfiletmp.close();
       remove(tmpNrgFilename.c_str());
       pfiletmp.setFileName(tmpNrgFilename);
     }
     explicit NrgData(nrgcore_type      *t_nrg_object, // nrgcore_type
                      const std::string &tfilename = "")
         : nrg_object(t_nrg_object) {
       if (tfilename.empty()) {
         tmpNrgFilename = "tempfile" + std::to_string(std::rand()) + ".h5";
       } else {
         tmpNrgFilename = tfilename;
       }
       pfiletmp.setFileName(tmpNrgFilename);
       // TODO(sp): clear the input operator
     }
     void setNRGObject(nrgcore_type *t_nrg_object) { nrg_object = t_nrg_object; }
     // ~NrgData() { // clear/delete the temp file
     //   clear();
     // }
     void saveFinalState() {
       pfiletmp.write<size_t>(savedNRGIndex, "savedNRGIndex");
       pfiletmp.write(nrg_object->relativeGroundStateEnergy,
                      "relativeGroundStateEnergy");
     }
     void close() {
       // save the nrg Index
       pfiletmp.close();
       isClosed = true;
     }
     void clear() {
       if (!isClosed) {
         this->close();
         remove(tmpNrgFilename.c_str());
         isClosed = true;
       }
     }
     void saveCurrentData() {
       // Things to save
       // current_hamiltonQ; // next hamiltonians
       // current_sysmQ;     // next symmetries
       // pre_sysmQ;         // previous symmetries
       // eigenvaluesQ;      // Eigenvalues
       // coupled_nQ_index;
       savedNRGIndex.push_back(nrg_object->nrg_iterations_cnt);
       std::string hgroup{"/NrgItr" +
                          std::to_string(nrg_object->nrg_iterations_cnt) + "/"};
       if (debugIO) {
         std::cout << "Writing file: " << tmpNrgFilename
                   << " for the Group:: " << hgroup << std::endl;
       }
       pfiletmp.createGroup(hgroup);
       // save the data
       pfiletmp.write<double>(nrg_object->current_hamiltonQ,
                              hgroup + "current_hamiltonQ");
       pfiletmp.write<int>(nrg_object->current_sysmQ, hgroup + "current_sysmQ");
       pfiletmp.write<int>(nrg_object->pre_sysmQ, hgroup + "pre_sysmQ");
       pfiletmp.write<double>(nrg_object->eigenvaluesQ, hgroup + "eigenvaluesQ");
       pfiletmp.write(nrg_object->coupled_nQ_index, hgroup + "coupled_nQ_index");
       pfiletmp.write<size_t>(nrg_object->eigenvaluesQ_kept_indices,
                              hgroup + "eigenvaluesQ_kept_indices");
       pfiletmp.write(nrg_object->all_eigenvalue,
                      "Eigenvalues" +
                          std::to_string(nrg_object->nrg_iterations_cnt));
       // Close the file
       //--------------------------------------------------------------
       // End of saveNrgData0
     }
     void saveqOperator(std::vector<qOperator> *opr, const std::string &hgroup) {
       // std::string hgroup{oprString};
       if (debugIO) {
         std::cout << "Writing file: " << tmpNrgFilename
                   << " for the Group:: " << hgroup << std::endl;
       }
       auto dg = pfiletmp.createGroup(hgroup);
       { dg.write_atr<size_t>(opr->size(), "operatorSize"); }
       // save the block size
       for (size_t ip = 0; ip < opr->size(); ip++) {
         auto localGr = hgroup + "/at" + std::to_string(ip) + "/";
         auto ds      = pfiletmp.createGroup(localGr);
         //
         std::vector<std::vector<size_t>> idVector;
         std::vector<size_t>              colVector;
         std::vector<size_t>              rowVector;
         size_t                           ic = 0;
         for (const auto &aa : *opr->at(ip).getMap()) {
           idVector.push_back({aa.first[0], aa.first[1]});
           colVector.push_back(aa.second.getcolumn());
           rowVector.push_back(aa.second.getrow());
           pfiletmp.write<double>(localGr + "qmat" + std::to_string(ic),
                                  aa.second.data(), aa.second.size());
           ic++;
         }
         pfiletmp.write(idVector, localGr + "idVector");
         pfiletmp.write<size_t>(colVector, localGr + "colVector");
         pfiletmp.write<size_t>(rowVector, localGr + "rowVector");
       }
       //--------------------------------------------------------------
       // End of saveNrgData0
     }
     void loadqOperator(std::vector<qOperator> *opr, const std::string &hgroup) {
       // clear the operator list
       opr->clear();
       //
       if (debugIO) {
         std::cout << "Reading file: " << tmpNrgFilename
                   << " for the Group:: " << hgroup << std::endl;
       }
       size_t operatorSize{0};
       {
         auto ds = pfiletmp.getGroup(hgroup);
         ds.read_atr(operatorSize, "operatorSize");
       }
       // save the block size
       opr->resize(operatorSize);
       for (size_t ip = 0; ip < operatorSize; ip++) {
         auto      localGr = hgroup + "/at" + std::to_string(ip) + "/";
         qOperator iOperator;
         //
         std::vector<std::vector<size_t>> idVector;
         std::vector<size_t>              colVector;
         std::vector<size_t>              rowVector;
         pfiletmp.read(idVector, localGr + "idVector");
         pfiletmp.read(colVector, localGr + "colVector");
         pfiletmp.read(rowVector, localGr + "rowVector");
         for (size_t ic = 0; ic < idVector.size(); ic++) {
           std::vector<double> qmat;
           pfiletmp.read<double, std::vector>(qmat, localGr + "qmat" +
                                                        std::to_string(ic));
           iOperator.set(qmatrix(qmat, rowVector[ic], colVector[ic]),
                         idVector[ic][0], idVector[ic][1]);
         }
         opr->at(ip) = iOperator;
       }
       //--------------------------------------------------------------
       // End of saveNrgData0
     }
     void loadCurrentData() { loadCurrentData(nrg_object->nrg_iterations_cnt); }
     void loadCurrentData(int in) {
       // Things to save
       // current_hamiltonQ; // next hamiltonians
       // current_sysmQ;     // next symmetries
       // pre_sysmQ;         // previous symmetries
       // eigenvaluesQ;      // Eigenvalues
       // coupled_nQ_index;
       nrg_object->nrg_iterations_cnt = in;
       std::string hgroup{"/NrgItr" + std::to_string(in) + "/"};
       if (debugIO) {
         std::cout << "Reading file: " << tmpNrgFilename
                   << " for the Group:: " << hgroup << std::endl;
       }
       // save the data
       pfiletmp.read(nrg_object->current_hamiltonQ, hgroup + "current_hamiltonQ");
       pfiletmp.read(nrg_object->current_sysmQ, hgroup + "current_sysmQ");
       pfiletmp.read(nrg_object->pre_sysmQ, hgroup + "pre_sysmQ");
       pfiletmp.read(nrg_object->eigenvaluesQ, hgroup + "eigenvaluesQ");
       pfiletmp.read(nrg_object->coupled_nQ_index, hgroup + "coupled_nQ_index");
       pfiletmp.read(nrg_object->eigenvaluesQ_kept_indices,
                     hgroup + "eigenvaluesQ_kept_indices");
       // if (save_f_operators) {
       //  loadqOperator(nrg_object->getWilsonSiteOperators(), "fdag_operator");
       //}
       //
       //--------------------------------------------------------------
       // End of loadNrgData0
     }
     bool                debugIO = false;
     std::vector<size_t> savedNRGIndex;
   };
