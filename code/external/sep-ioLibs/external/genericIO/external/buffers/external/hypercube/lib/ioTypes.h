#ifndef iotypes_h
#define iotypes_h 1
#include <string>
namespace SEP {

/*!
enum for specifyting the type of data*/
enum dataType {
  DATA_BYTE,     /// 8-bit unsigned char
  DATA_INT,      /// 4-byte integer
  DATA_FLOAT,    /// 4-byte float
  DATA_COMPLEX,  /// 2 4-byte floats using c++11 templates
  DATA_DOUBLE,   /// 8-byte float
  DATA_SHORT,    /// 2-byte integer
  DATA_UNKNOWN   /// Unspecifed type
};

size_t getDataTypeSize(
    const dataType typ);  /// Return the size of the given type
dataType toElementType(const std::string &name);  /// From string return
                                                  /// dataType
std::string getTypeString(const dataType typ);  /// From data type return na,e

}  // namespace SEP
#endif
