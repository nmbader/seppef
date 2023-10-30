

#include "ioTypes.h"
#include <complex>
#include <string>
#include "SEPException.h"
using namespace SEP;
//! Convert string to element type
/*!
\param name String name of element type */
dataType SEP::toElementType(const std::string &name) {
  if (name == "BYTE") return DATA_BYTE;
  if (name == "INT") return DATA_INT;
  if (name == "FLOAT") return DATA_FLOAT;
  if (name == "COMPLEX") return DATA_COMPLEX;
  if (name == "COMPLEXDOUBLE") return DATA_COMPLEXDOUBLE;
  if (name == "DOUBLE") return DATA_DOUBLE;
  if (name == "SHORT") return DATA_SHORT;
  return DATA_UNKNOWN;
}
//! Convert element type to string
/*!
\param typ data element type */
std::string SEP::getTypeString(const dataType typ) {
  switch (typ) {
    case DATA_BYTE:
      return "BYTE";
      break;
    case DATA_FLOAT:
      return "FLOAT";
      break;
    case DATA_COMPLEX:
      return "COMPLEX";
      break;
    case DATA_INT:
      return "INT";
      break;
    case DATA_SHORT:
      return "SHORT";
      break;
    case DATA_DOUBLE:
      return "DOUBLE";
      break;
    case DATA_COMPLEXDOUBLE:
      return "COMPLEXDOUBLE";
      break;
    default:
      throw(SEPException(std::string("Unknown data type")));
  }
}
//! Get the size of element type
/*!
\param typ Data element type */
size_t SEP::getDataTypeSize(const dataType typ) {
  switch (typ) {
    case DATA_BYTE:
      return sizeof(unsigned char);
      break;
    case DATA_FLOAT:
      return sizeof(float);
      break;
    case DATA_COMPLEXDOUBLE:
      return sizeof(std::complex<double>);
      break;
    case DATA_COMPLEX:
      return sizeof(std::complex<float>);
      break;
    case DATA_INT:
      return sizeof(int);
      break;
    case DATA_DOUBLE:
      return sizeof(double);
      break;
    default:
      throw(SEPException(std::string("Unknown data type")));
  }
}
