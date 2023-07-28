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
/**
 * @brief This class provides a way to save the  `nrgcore` internal  data
 * from the NRG calculation
 * into a  HDF5 file. The data is saved in a file and can be read back
 * into the NRG calculation. This is useful for back-ward iteration of the
 * of the nrg-iterations. Everything needed for NRG iteration can accessed
 * from this this file. This class also provides functionality to
 * read and write `qOperator`.
 *
 *
 *
 * @tparam nrgcore_type Type of the `nrgcore` object.
 * @param tfilename File name of the HDF5 file to save the data. If this
 * is empty, a random file name is generated starting with the name  `tempfile`.
 */
template <typename nrgcore_type> class NrgData {
  std::string        tmpNrgFilename;
  nrgcore_type      *nrg_object;
  h5stream::h5stream pfiletmp;
  bool               isClosed = false;

public:
  /**
   * @brief This sets the `nrg_object` to `nullptr` and you need to specify
   *  the `nrg_object` by calling the `setNRGObject` function.
   *
   * @param tfilename File name of the HDF5 file to save the data.
   */
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
  /**
   * @brief Read the data from the file `tfile` into the `nrg_object`. This
   * function is used to read the data back into the `nrg_object` for back-ward
   * iteration of the NRG calculation.
   *
   *
   */
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
  /**
   * @brief Set the file name of the HDF5 file to save the data.
   *
   */
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
  /**
   * @brief Set the `nrg_object` to `t_nrg_object`.
   *
   */
  void setNRGObject(nrgcore_type *t_nrg_object) { nrg_object = t_nrg_object; }
  // ~NrgData() { // clear/delete the temp file
  //   clear();
  // }
  /**
   * @brief We should call this function after the last iteration of the Wilson
   * Chain.
   */
  void saveFinalState() {
    pfiletmp.write<size_t>(savedNRGIndex, "savedNRGIndex");
    pfiletmp.write(nrg_object->relativeGroundStateEnergy,
                   "relativeGroundStateEnergy");
  }
  /**
   * @brief We should always call this function after the last iteration of the
   * Wilson Chain to close file.
   */
  void close() {
    // save the nrg Index
    pfiletmp.close();
    isClosed = true;
  }
  /**
   * @brief This \f \textbf{deletes} \f the temp file and closes the file.
   */
  void clear() {
    if (!isClosed) {
      this->close();
      remove(tmpNrgFilename.c_str());
      isClosed = true;
    }
  }
  /**
   * @brief This saves data for the  the current iteration state
   * of the Wilson Chain.
   */
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
  /**
   * @brief Saves the `qOperator` in the file. This function is called
   * after the `qOperator` is constructed and rotated in the eigenbasis
   * after each iteration.
   *
   * @param opr  `qOperator` to be saved.
   * @param hgroup The group name in which the `qOperator` is to be saved.
   * Remember pass same string when reading the `qOperator` from the file.
   */
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
  /**
   * @brief  Reads the `qOperator` from the file. This function is called
   * in the back-ward iteration of the Wilson Chain.
   *
   *
   * @param opr  `qOperator` to be read from the file.
   * @param hgroup The group name in which the `qOperator` was saved.
   */
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
  /**
   * @brief Loads the data from the file. This function is called in the
   * back-ward iteration of the Wilson Chain.
   *
   * @param in The current iteration number or the Wilson number. We normally
   * pass the `nrg_object->nrg_iterations_cnt` to this function. We set the
   * impurity to the iteration number to be -1.
   */
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
  /**
   * @brief Print Debugg Info
   */
  bool debugIO = false;
  /**
   * @brief Nrg iteration Numbers for which higher energy states are discarded.
   */
  std::vector<size_t> savedNRGIndex;
};
