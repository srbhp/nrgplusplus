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
 * @class NrgData
 * @brief Persist NRG iteration state and operators in an HDF5 file.
 *
 * This class stores the snapshot of a Wilson-chain iteration needed for later
 * backward iteration, including the current Hamiltonians, symmetry sectors,
 * eigenvalues, kept indices, and any `qOperator` blocks. The data can be saved
 * to disk, reloaded, and reused when reconstructing the chain state.
 *
 * @tparam nrgcore_type Type of the underlying NRG core object.
 */
template <typename nrgcore_type> class NrgData {
  std::string        tmpNrgFilename;
  nrgcore_type      *nrg_object;
  h5stream::h5stream pfiletmp;
  bool               isClosed = false;

public:
  /**
   * @brief Construct a data container without attaching an NRG object yet.
   *
   * The object may be associated with a concrete NRG core later using
   * `setNRGObject`.
   *
   * @param tfilename Optional path to the HDF5 file to be used for persistence.
   * If empty, a temporary file name is generated automatically.
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
   * @brief Read previously saved iteration data from disk into the active NRG
   * object.
   *
   * This is typically used during backward iteration to restore a stored shell
   * state from the HDF5 archive.
   *
   * @param tfile Path to the HDF5 file containing the saved data.
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
   * @brief Update the HDF5 file name used for persistence.
   *
   * @param tfile New output file path.
   */
  void setFileName(const std::string &tfile) {
    tmpNrgFilename = tfile;
    pfiletmp.close();
    remove(tmpNrgFilename.c_str());
    pfiletmp.setFileName(tmpNrgFilename);
  }

  /**
   * @brief Construct a data container bound to a specific NRG core object.
   *
   * @param t_nrg_object Pointer to the active NRG model state.
   * @param tfilename Optional HDF5 file path. If empty, a temporary file is
   * created automatically.
   */
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
   * @brief Attach an NRG core instance to this data container.
   *
   * @param t_nrg_object Pointer to the NRG object whose state will be stored and
   * restored.
   */
  void setNRGObject(nrgcore_type *t_nrg_object) { nrg_object = t_nrg_object; }

  /**
   * @brief Save the final accumulated state metadata after the last iteration.
   *
   * This stores the iteration index list and the final relative ground-state
   * energy for later reconstruction steps.
   */
  void saveFinalState() {
    pfiletmp.write<size_t>(savedNRGIndex, "savedNRGIndex");
    pfiletmp.write(nrg_object->relativeGroundStateEnergy,
                   "relativeGroundStateEnergy");
  }

  /**
   * @brief Close the HDF5 file associated with this object.
   */
  void close() {
    // save the nrg Index
    pfiletmp.close();
    isClosed = true;
  }

  /**
   * @brief Delete the temporary HDF5 file and close the stream.
   *
   * This removes the persistent snapshot associated with the current run.
   */
  void clear() {
    if (!isClosed) {
      this->close();
      remove(tmpNrgFilename.c_str());
      isClosed = true;
    }
  }

  /**
   * @brief Save the full state of the current Wilson-chain iteration.
   *
   * The function serializes the current shell Hamiltonians, symmetry sectors,
   * eigenvalues, coupling metadata, and kept-state indices into the output HDF5
   * file under a unique iteration group.
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
   * @brief Save a `qOperator` block into the HDF5 archive.
   *
   * @param opr Pointer to the operator list to write.
   * @param hgroup Name of the HDF5 group under which the operator data is stored.
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
   * @brief Restore a `qOperator` block from the HDF5 archive.
   *
   * @param opr Pointer to the operator list to populate.
   * @param hgroup Name of the HDF5 group containing the saved operator data.
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

  /**
   * @brief Load the currently active iteration from the file.
   *
   * This convenience wrapper calls `loadCurrentData(nrg_object->nrg_iterations_cnt)`.
   */
  void loadCurrentData() { loadCurrentData(nrg_object->nrg_iterations_cnt); }

  /**
   * @brief Restore a saved shell state from the HDF5 archive.
   *
   * @param in Iteration index to load from the stored data set.
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
   * @brief Enable or disable verbose debug logging for HDF5 I/O operations.
   */
  bool debugIO = false;

  /**
   * @brief List of saved NRG iteration indices retained in the archive.
   */
  std::vector<size_t> savedNRGIndex;
};
