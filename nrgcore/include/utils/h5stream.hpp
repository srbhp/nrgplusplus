/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once
#include "H5Cpp.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
// #include <vector>
//  Maps C++ type to HDF5 type
template <typename T> inline const H5::PredType &get_datatype_for_hdf5();
// Reference:
// https://www.hdfgroup.org/HDF5/doc/cpplus_RM/class_h5_1_1_pred_type.html
template <> inline const H5::PredType &get_datatype_for_hdf5<char>() {
  return H5::PredType::NATIVE_CHAR;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned char>() {
  return H5::PredType::NATIVE_UCHAR;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<short>() {
  return H5::PredType::NATIVE_SHORT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned short>() {
  return H5::PredType::NATIVE_USHORT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<int>() {
  return H5::PredType::NATIVE_INT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned int>() {
  return H5::PredType::NATIVE_UINT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<long>() {
  return H5::PredType::NATIVE_LONG;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned long>() {
  return H5::PredType::NATIVE_ULONG;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<long long>() {
  return H5::PredType::NATIVE_LLONG;
}
template <>
inline const H5::PredType &get_datatype_for_hdf5<unsigned long long>() {
  return H5::PredType::NATIVE_ULLONG;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<float>() {
  return H5::PredType::NATIVE_FLOAT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<double>() {
  return H5::PredType::NATIVE_DOUBLE;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<long double>() {
  return H5::PredType::NATIVE_LDOUBLE;
}
//---------------------------------------------------------
namespace h5stream {
template <typename T> struct h5str1 {
  std::string    keyName;
  T             *data;
  const unsigned dataSize{};
};
} // namespace h5stream
//
namespace h5stream {
/**
 * @class dspace
 * @brief This class is responsible to read and write of the metadata/atrributes
 * fo a dataset.
 *
 */
class dspace {
public:
  H5::DataSet dataset;
  explicit dspace(const H5::DataSet &datasetx) : dataset(datasetx){};
  /**
   * @brief Write Metadata
   *
   * @tparam T
   * @param data
   * @param dataname
   */
  template <typename T>
  void write_atr(const T data, const H5std_string &dataname) {
    auto          type           = get_datatype_for_hdf5<T>();
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attribute =
        dataset.createAttribute(dataname, type, attr_dataspace);
    attribute.write(type, &data);
  }
  /**
   * @brief Read Metadata
   *
   * @tparam T
   * @param data
   * @param dataname
   */
  template <typename T> void read_atr(T &data, const H5std_string &dataname) {
    // auto type = get_datatype_for_hdf5<T>();
    H5::Attribute attribute = dataset.openAttribute(dataname);
    H5::DataType  type      = attribute.getDataType();
    attribute.read(type, &data);
  }
  //----------------------
};
// ---- -----------------------------------
/**
 * @class gspace: Read and write metadata to  a `H5::Group` object.
 * @brief
 *
 */
class gspace { // for group
public:
  H5::Group dataset;
  explicit gspace(const H5::Group &datasetx) : dataset(datasetx){};
  /**
   * @brief Write Metadata
   *
   * @tparam T
   * @param data
   * @param dataname
   */
  template <typename T>
  void write_atr(const T data, const H5std_string &dataname) {
    auto          type           = get_datatype_for_hdf5<T>();
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attribute =
        dataset.createAttribute(dataname, type, attr_dataspace);
    attribute.write(type, &data);
  }
  /**
   * @brief Read Metadata
   *
   * @tparam T
   * @param data
   * @param dataname
   */
  template <typename T> void read_atr(T &data, const H5std_string &dataname) {
    // auto type = get_datatype_for_hdf5<T>();
    H5::Attribute attribute = dataset.openAttribute(dataname);
    H5::DataType  type      = attribute.getDataType();
    attribute.read(type, &data);
  }
  //----------------------
};
} // namespace h5stream
//******************************************************************
// namespace h5stream
//----------------------------------------------------------
namespace h5stream {
/**
 * @class h5stream::h5stream
 * @brief C++ Header-only library for simple HDF5 input/output
 *
 *  \section <1>  How to use
 *
 *		Just include the `h5stream.hpp` into your your main file.
 *
 *  \subsection <1.1> How to  Compile
 *
 *		@code{.bash}
 *		g++ -lhdf5 -lhdf5_cpp -std=c++1z example.cpp
 *		@endcode
 *
 * \section <2.1> Example
 *
 *  \subsection <2.2> Create a File with a mode.
 *
 *
 *		- "tr":   Create a file, truncate if it exists, Default
 *		- "r":    Readonly, the file must exist
 *		- "rw": Read/write, the file must exist
 *		- "x":   Create a file, fail if exists
 *
 *
 *		@code{.cpp}
 *		h5stream::h5stream file("sample.h5", "tr");
 *		// or
 *		h5stream::h5stream file("sample.h5");
 *		@endcode
 *
 * \subsection <2.3> write and read `std::vector`
 *
 *		Create a vector and write it to the file
 *
 *
 *		@code{.cpp}
 *		std::vector<double> matrix { 1, 2, 3282, 932 };
 *		file.write<double>(matrix, "matrix");
 *		@endcode
 *
 *
 * \subsection <2.4>  write and read Metadata
 *
 *		Write Attributes( Metadata) to the to the same data space
 *
 *
 *		@code{.cpp}
 *		auto dspace = file.get_dataspace("matrix");
 *		dspace.write_atr<double>(1.2, "Units");
 *		@endcode
 *
 *
 * \subsection <2.5> Read data from the file
 *
 *
 *		@code{.cpp}
 *		auto xx = file.read_vector<double>("matrix");
 *		//OR
 *		file.read<double>(xx, "matrix");
 *		@endcode
 *
 *
 *  \subsection <2.6>  Read Attribute (Metadata)
 *		@code{.cpp}
 *		double x = 0;
 *		dspace.read_atr<double>(x, "Units");
 *		std::cout << "Attribute : " << x << std::endl;
 *		std::cout << "HDF file size (MB): " << file.file_size() <<
 *       std::endl;
 *		@endcode
 *
 */
class h5stream {
  bool debug = false;

public:
  H5std_string hdf5FileName;
  H5::H5File   hdf5File;
  h5stream() = default;
  ;
  explicit h5stream(const std::string &fileName,
                    const std::string &rw = std::string("tr")) {
    setFileName(fileName, rw);
  }
  /**
   * @brief create and set the file name
   *
   * @param fileName
   * @param rw : Possible values are
   * "r": read only,
   *
   * "rw": Read-write access. If the file is currently open for
   * read-only access then it will be reopened. Absence of this
   * flag implies read-only access.
   *
   * "x": Fail if file already exists. "r" and "tr" are mutually exclusive
   *
   * "tr": Truncate file, if it already exists, erasing all data previously
   * stored in the file.: default
   */
  void setFileName(const H5std_string &fileName, // NOLINT
                   const std::string  &rw = std::string("tr")) {
    hdf5FileName = fileName;
    std::cout << "hdf5FileName:" << hdf5FileName << std::endl;
    try {
      H5::Exception::dontPrint();
      if (rw == "r") {
        hdf5File = H5::H5File(fileName, H5F_ACC_RDONLY);
      }
      if (rw == "rw") {
        hdf5File = H5::H5File(fileName, H5F_ACC_RDWR);
      }
      if (rw == "x") {
        hdf5File = H5::H5File(fileName, H5F_ACC_EXCL);
      }
      if (rw == "tr") {
        hdf5File = H5::H5File(fileName, H5F_ACC_TRUNC);
      }
    } catch (...) {
      std::cout << " Error :: Unable to setFileName!  " << fileName
                << std::endl;
    }
  }
  template <typename T = double, template <typename...> class vec>
  void write(const std::vector<vec<T>> &data, const H5std_string &datasetName) {
    for (size_t i = 0; i < data.size(); i++) {
      write<T, vec>(data[i], datasetName + std::to_string(i));
    }
  }
  /**
   * @brief
   *
   * @param data
   * @param datasetName
   */
  template <typename T = double, template <typename...> class vec = std::vector>
  void write(const vec<T> &data, const H5std_string &datasetName) {
    write<T>(datasetName, data.data(), data.size());
  }
  // Write raw pointer
  /**
   * @brief
   *
   * @param datasetName
   * @param data
   * @param data_size
   */
  template <typename T = double>
  void write(const H5std_string &datasetName, const T *data,
             unsigned data_size) {
    try {
      H5::Exception::dontPrint();
      const int RANK = 1;
      auto      type = get_datatype_for_hdf5<T>();
      hsize_t   dimsf[1];   // dataset dimensions
      dimsf[0] = data_size; //
      H5::DataSpace dataspace(RANK, dimsf);
      H5::DataSet   dataset =
          hdf5File.createDataSet(datasetName, type, dataspace);
      if (data_size != 0) { // Dont write if it zero
        dataset.write(data, type);
      }
    } catch (...) {
      std::string errString = "Error! :: Unable to write datasetName " +
                              datasetName + " in  the file " + hdf5FileName;
      throw std::runtime_error(errString);
    }
  }
  // Read file
  /**
   * @brief
   *
   * @param data
   * @param datasetName
   */
  template <typename T = double, template <typename...> class vec = std::vector>
  void read(vec<T> &data, const H5std_string &datasetName) {
    try {
      H5::Exception::dontPrint();
      auto          type      = get_datatype_for_hdf5<T>();
      H5::DataSet   dataset   = hdf5File.openDataSet(datasetName);
      H5::DataSpace dataspace = dataset.getSpace();
      hsize_t       dim[1];
      dataspace.getSimpleExtentDims(dim, nullptr);
      data.resize(dim[0]);
      if (dim[0] != 0) { // Dont read if it zero
        dataset.read(data.data(), type, dataspace, dataspace);
      }
    } catch (...) {
      std::string errString = "Error! :: Unable to READ datasetName " +
                              datasetName + " from the file " + hdf5FileName;
      throw std::runtime_error(errString);
      //    H5::FileIException::printErrorStack();
    }
  }
  // Higher order vector  or matrix type
  // User has to give the correct size of the std::vector
  /**
   * @brief Read the data
   *
   * @param data
   * @param datasetName
   */
  template <typename T = double, template <typename...> class vec>
  void read(std::vector<vec<T>> &data, const H5std_string &datasetName) {
    data.clear();
    size_t icount{0};
    bool   foundKey{true};
    while (foundKey) {
      // for (size_t i = 0; i < 10; i++) {
      try {
        std::vector<T> aa;
        read<T, std::vector>(aa, datasetName + std::to_string(icount));
        data.push_back(aa);
        icount++;
      } catch (std::exception &e) {
        // foundKey = false;
        // std::cout << "foundKey: " << foundKey << e.what() << std::endl;
        break;
      }
    }
    // throw warning if the dataset name was not found
    if (icount == 0) {
      std::string err_string =
          "Error :: Unable to read dataset! " + std::string(datasetName);
      throw std::runtime_error(err_string);
    }
  }
  /**
   * @brief close the file
   */
  void close() { hdf5File.close(); }
  /**
   * @brief Return the filesize in MB.
   *
   * @return
   */
  [[nodiscard]] double fileSize() const {
    return static_cast<double>(hdf5File.getFileSize()) / (1024 * 1024.);
  }
  /**
   * @brief
   *
   * @param dataset_name
   * @return
   */
  dspace getDataspace(const H5std_string &dataset_name) {
    return dspace(hdf5File.openDataSet(dataset_name));
  }
  /**
   * @brief
   *
   * @param dataset_name
   * @return
   */
  gspace getGroup(const H5std_string &dataset_name) {
    return gspace(hdf5File.openGroup(dataset_name));
  }
  /**
   * @brief Create a group name
   *
   * @param group_name
   * @return
   */
  auto createGroup(const H5std_string &group_name) {
    return gspace(hdf5File.createGroup(group_name));
  }
  //*************************************************
  // write metadata at the root label
  /**
   * @brief Write Metadata
   *
   * @tparam T `datatype`
   * @param data
   * @param label
   */
  template <typename T>
  void writeMetadata(const T &data, const H5std_string &label) {
    //
    auto ds = getDataspace("");
    ds.write_atr(data, label);
  }
  /**
   * @brief
   *
   * @tparam T
   * @param data
   * @param label
   */
  template <typename T> void readMetadata(T &data, const H5std_string &label) {
    //
    auto ds = getDataspace("");
    ds.read_atr(data, label);
  }
  //*************************************************
  // overload << and >>
  template <typename T>
  friend h5stream &operator<<(h5stream &out, const h5str1<T> &struct1) {
    out.write<T>(struct1.keyName, struct1.data, struct1.dataSize);
    return out;
  }
  template <typename T>
  friend h5stream &operator>>(h5stream &out, const h5str1<T> &struct1) {
    out.read<T>(struct1.keyName, struct1.data, struct1.dataSize);
    return out;
  }
  //-------------------------------------------------
};
} // namespace h5stream
