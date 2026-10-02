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
#include <H5Cpp.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Maps a C++ type to the corresponding HDF5 data type.
 *
 * This function template is specialized for various C++ types to return
 * the appropriate HDF5 data type.
 *
 * @tparam T The C++ type to map.
 * @return A reference to the corresponding HDF5 data type.
 */
template <typename T>
inline const H5::PredType &get_datatype_for_hdf5();

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
 * @brief Helper for reading and writing dataset attributes in HDF5 files.
 *
 * This lightweight wrapper exposes attribute read/write operations bound to an
 * existing `H5::DataSet` instance.
 */
class dspace {
public:
  H5::DataSet dataset;

  /**
   * @brief Construct a dataset wrapper around an existing HDF5 dataset.
   *
   * @param datasetx Underlying HDF5 dataset object.
   */
  explicit dspace(const H5::DataSet &datasetx) : dataset(datasetx) {};

  /**
   * @brief Write a scalar attribute to the wrapped dataset.
   *
   * @tparam T Type of the metadata value.
   * @param data Scalar value to store as an attribute.
   * @param dataname Name of the attribute.
   */
  template <typename T>
  void write_atr(const T &data, const H5std_string &dataname) {
    auto          type           = get_datatype_for_hdf5<T>();
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attribute =
        dataset.createAttribute(dataname, type, attr_dataspace);
    attribute.write(type, &data);
  }

  /**
   * @brief Read a scalar attribute from the wrapped dataset.
   *
   * @tparam T Type of the metadata value.
   * @param data Value to populate from the attribute.
   * @param dataname Name of the attribute.
   */
  template <typename T> void read_atr(T &data, const H5std_string &dataname) {
    H5::Attribute attribute = dataset.openAttribute(dataname);
    H5::DataType  type      = attribute.getDataType();
    attribute.read(type, &data);
  }
};

/**
 * @class gspace
 * @brief Helper for reading and writing attributes on an HDF5 group.
 *
 * This wrapper mirrors `dspace` but binds directly to an HDF5 group instead of
 * a dataset.
 */
class gspace {
public:
  H5::Group dataset;

  /**
   * @brief Construct a group wrapper around an existing HDF5 group.
   *
   * @param datasetx Underlying HDF5 group object.
   */
  explicit gspace(const H5::Group &datasetx) : dataset(datasetx) {};

  /**
   * @brief Write a scalar attribute to the wrapped group.
   *
   * @tparam T Type of the metadata value.
   * @param data Scalar value to store as an attribute.
   * @param dataname Name of the attribute.
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
   * @brief Read a scalar attribute from the wrapped group.
   *
   * @tparam T Type of the metadata value.
   * @param data Value to populate from the attribute.
   * @param dataname Name of the attribute.
   */
  template <typename T> void read_atr(T &data, const H5std_string &dataname) {
    H5::Attribute attribute = dataset.openAttribute(dataname);
    H5::DataType  type      = attribute.getDataType();
    attribute.read(type, &data);
  }
};
} // namespace h5stream

namespace h5stream {
/**
 * @class h5stream
 * @brief Minimal header-only wrapper for common HDF5 read/write operations.
 *
 * The class provides utilities to create or open an HDF5 file, write/read
 * vectors and scalar data, attach metadata attributes, and navigate groups.
 */
class h5stream {
  bool debug = false;

public:
  H5std_string hdf5FileName;
  H5::H5File   hdf5File;

  /**
   * @brief Default construct an empty HDF5 stream object.
   */
  h5stream() = default;

  /**
   * @brief Construct and immediately open an HDF5 file.
   *
   * @param fileName Path to the HDF5 file.
   * @param rw Access mode. Supported values are `"tr"`, `"r"`, `"rw"`, and
   * `"x"`.
   */
  explicit h5stream(const std::string &fileName,
                    const std::string &rw = std::string("tr")) {
    setFileName(fileName, rw);
  }

  /**
   * @brief Open or recreate the HDF5 file with the requested mode.
   *
   * @param fileName Path to the file.
   * @param rw Access mode: `"r"`, `"rw"`, `"x"`, or `"tr"`.
   */
  void setFileName(const H5std_string &fileName,
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

  /**
   * @brief Write a nested vector structure to sequential HDF5 datasets.
   *
   * @tparam T Element type of the vector.
   * @tparam vec Container type used for the nested vector.
   * @param data Vector of vectors to serialize.
   * @param datasetName Base name for the generated datasets.
   */
  template <typename T = double, template <typename...> class vec>
  void write(const std::vector<vec<T>> &data, const H5std_string &datasetName) {
    for (size_t i = 0; i < data.size(); i++) {
      write<T, vec>(data[i], datasetName + std::to_string(i));
    }
  }

  /**
   * @brief Write a single vector to an HDF5 dataset.
   *
   * @tparam T Element type.
   * @tparam vec Container type used for the vector.
   * @param data Vector payload.
   * @param datasetName Dataset name in the file.
   */
  template <typename T = double, template <typename...> class vec = std::vector>
  void write(const vec<T> &data, const H5std_string &datasetName) {
    write<T>(datasetName, data.data(), data.size());
  }

  /**
   * @brief Write a raw C-style array to an HDF5 dataset.
   *
   * @tparam T Element type.
   * @param datasetName Dataset name in the HDF5 file.
   * @param data Pointer to the first element.
   * @param data_size Number of elements in the array.
   */
  template <typename T = double>
  void write(const H5std_string &datasetName, const T *data,
             unsigned data_size) {
    try {
      H5::Exception::dontPrint();
      const int RANK = 1;
      auto      type = get_datatype_for_hdf5<T>();
      hsize_t   dimsf[1];
      dimsf[0] = data_size;
      H5::DataSpace dataspace(RANK, dimsf);
      H5::DataSet   dataset =
          hdf5File.createDataSet(datasetName, type, dataspace);
      if (data_size != 0) {
        dataset.write(data, type);
      }
    } catch (...) {
      std::string errString = "Error! :: Unable to write datasetName " +
                              datasetName + " in  the file " + hdf5FileName;
      throw std::runtime_error(errString);
    }
  }

  /**
   * @brief Read a vector from an HDF5 dataset.
   *
   * @tparam T Element type.
   * @tparam vec Container type used for the vector.
   * @param data Output container to populate.
   * @param datasetName Name of the dataset in the file.
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
      if (dim[0] != 0) {
        dataset.read(data.data(), type, dataspace, dataspace);
      }
    } catch (...) {
      std::string errString = "Error! :: Unable to READ datasetName " +
                              datasetName + " from the file " + hdf5FileName;
      throw std::runtime_error(errString);
    }
  }

  /**
   * @brief Read a nested vector structure stored as sequential datasets.
   *
   * @tparam T Element type.
   * @tparam vec Container type used for the nested vector.
   * @param data Output vector of vectors.
   * @param datasetName Base name for the generated datasets.
   */
  template <typename T = double, template <typename...> class vec>
  void read(std::vector<vec<T>> &data, const H5std_string &datasetName) {
    data.clear();
    size_t icount{0};
    bool   foundKey{true};
    while (foundKey) {
      try {
        std::vector<T> aa;
        read<T, std::vector>(aa, datasetName + std::to_string(icount));
        data.push_back(aa);
        icount++;
      } catch (std::exception &e) {
        break;
      }
    }
    if (icount == 0) {
      std::string err_string =
          "Error :: Unable to read dataset! " + std::string(datasetName);
      throw std::runtime_error(err_string);
    }
  }

  /**
   * @brief Close the currently open HDF5 file.
   */
  void close() { hdf5File.close(); }

  /**
   * @brief Return the file size in megabytes.
   *
   * @return File size in MB.
   */
  [[nodiscard]] double fileSize() const {
    return static_cast<double>(hdf5File.getFileSize()) / (1024 * 1024.);
  }

  /**
   * @brief Open and wrap a dataset for attribute access.
   *
   * @param dataset_name Name of the dataset.
   * @return Lightweight dataset attribute wrapper.
   */
  dspace getDataspace(const H5std_string &dataset_name) {
    return dspace(hdf5File.openDataSet(dataset_name));
  }

  /**
   * @brief Open and wrap a group for attribute access.
   *
   * @param dataset_name Name of the group.
   * @return Lightweight group attribute wrapper.
   */
  gspace getGroup(const H5std_string &dataset_name) {
    return gspace(hdf5File.openGroup(dataset_name));
  }

  /**
   * @brief Create a new HDF5 group.
   *
   * @param group_name Name of the group to create.
   * @return Wrapper around the newly created group.
   */
  auto createGroup(const H5std_string &group_name) {
    return gspace(hdf5File.createGroup(group_name));
  }

  /**
   * @brief Write a scalar metadata attribute at the root level.
   *
   * @tparam T Type of the metadata value.
   * @param data Metadata value to store.
   * @param label Attribute name.
   */
  template <typename T>
  void writeMetadata(const T &data, const H5std_string &label) {
    auto ds = getDataspace("");
    ds.write_atr(data, label);
  }

  /**
   * @brief Read a scalar metadata attribute from the root level.
   *
   * @tparam T Type of the metadata value.
   * @param data Output value to populate.
   * @param label Attribute name.
   */
  template <typename T> void readMetadata(T &data, const H5std_string &label) {
    auto ds = getDataspace("");
    ds.read_atr(data, label);
  }

  /**
   * @brief Stream-style write helper for a simple HDF5 record descriptor.
   *
   * @tparam T Type of the pointed data.
   * @param out Output stream wrapper.
   * @param struct1 Descriptor containing the dataset name and pointer.
   * @return Reference to the modified stream.
   */
  template <typename T>
  friend h5stream &operator<<(h5stream &out, const h5str1<T> &struct1) {
    out.write<T>(struct1.keyName, struct1.data, struct1.dataSize);
    return out;
  }

  /**
   * @brief Stream-style read helper for a simple HDF5 record descriptor.
   *
   * @tparam T Type of the pointed data.
   * @param out Input stream wrapper.
   * @param struct1 Descriptor containing the dataset name and pointer.
   * @return Reference to the modified stream.
   */
  template <typename T>
  friend h5stream &operator>>(h5stream &out, const h5str1<T> &struct1) {
    out.read<T>(struct1.keyName, struct1.data, struct1.dataSize);
    return out;
  }
};
} // namespace h5stream
