#ifndef HEADER_IO_H
#define HEADER_IO_H 1
#include <map>
#include <memory>
#include <vector>
#include "SEPException.h"
#include "ioTypes.h"

namespace SEP {
/*!
Basic key class. Just contains name and data type */
class key {
 public:
  //! Create key
  /*!
    \param name the name of the key
    \param typ the type of the key
  */
  key(const std::string name, const dataType typ) {
    _name = name;
    _type = typ;
  }
  //! Return dataType
  dataType type() const { return _type; }
  //! Return name of key

  std::string name() const { return _name; }

 private:
  std::string _name;  /// Name of the key
  dataType _type;     /// Datatype for key
};

/*!
 Class containing a collection of keys and values */
class header {
 public:
  //! Create header
  /*!
    \param keys List of keys associated with header
    \param nh Number of headers
  */
  header(const std::vector<key> keys, const int nh);

  //! Add a key to the headers
  /*!
    \param name  Name of of the key to add
    \param typ Type of key to add
  */
  void addKey(const std::string name, const dataType typ);
  //! Set all the values of a float key
  /*!
    \param name  Name of key
    \param vals Values to set
  */

  void setFloatKey(const std::string name, const std::vector<float> vals);
  //! Set all the values of a int key
  /*!
    \param name  Name of key
    \param vals Values to set
  */
  void setIntKey(const std::string name, const std::vector<int> vals);
  //! Set all the values of a double key
  /*!
    \param name  Name of key
    \param vals Values to set
  */
  void setDoubleKey(const std::string name, const std::vector<double> vals);
  //! Set all the values of a short key
  /*!
    \param name  Name of key
    \param vals Values to set
  */
  void setShortKey(const std::string name, const std::vector<short> vals);

  //! Copy in header values
  /*!
    \param head  All of the header values
  */
  void setHeaders(std::vector<std::vector<unsigned char>>& head);
  //! Get all float keys
  /*!
    \param name Name of the key to set
  */
  std::vector<float> getFloatKey(const std::string name) const;
  //! Get all double keys
  /*!
    \param name Name of the key to set
  */
  std::vector<double> getDoubleKey(const std::string name) const;
  //! Get all int keys
  /*!
    \param name Name of the key to set
  */
  std::vector<int> getIntKey(const std::string name) const;
  //! Get all short keys
  /*!
    \param name Name of the key to set
  */
  std::vector<short> getShortKey(const std::string name) const;
  //! Get an individual float key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
  */
  float getFloatKeyVal(const std::string name, const int index) const;
  //! Get an individual double key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
  */
  double getDoubleKeyVal(const std::string name, const int index) const;
  //! Get an individual int key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
  */
  int getIntKeyVal(const std::string name, const int index) const;
  //! Get an individual short key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
  */
  short getShortKeyVal(const std::string name, const int index) const;
  //! Set an individual float key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
        \param val Value to set key

  */
  void setFloatKeyVal(const std::string name, const int index, const float val);
  //! Set an individual short key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
    \param val Value to set key
  */
  void setShortKeyVal(const std::string name, const int index, const short val);
  //! Set an individual double key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
    \param val Value to set key
  */
  void setDoubleKeyVal(const std::string name, const int index,
                       const double val);
  //! Set an individual integer key
  /*!
    \param name Name of the key to set
    \param index Index of key to set
    \param val Value to set key
  */
  void setIntKeyVal(const std::string name, const int index, const int val);

  //! Clone headers
  std::shared_ptr<header> clone();

 private:
  //!  Get index of a given key
  /*!
    \param name Name of key to grab index of
  */
  size_t getIndex(std::string const name) const { return _key_index.at(name); }
  //! get offset of a given key
  /*!
    \param name Name of key to grab offset of
  */
  inline size_t getOffset(std::string const name) const {
    if (_key_offset.count(name) == 0)
      throw SEPException(std::string("Unknown key name") + name);

    return _key_offset.at(name);
  }

  std::vector<key> _keys;                     /// List of keys
  std::map<std::string, size_t> _key_index;   /// Index of each key
  std::map<std::string, size_t> _key_offset;  /// Offset (in byte) of givern key
  std::vector<std::vector<unsigned char>> _head;  /// Headers

  int _nsz;  /// Number of headers
};

}  // namespace SEP
#endif
