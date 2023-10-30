#ifndef file_buffers_reg_file_h
#define file_buffers_reg_file_h 1
#include <stdbool.h>
#include <string>
#include "basicIO.h"
#include "buffers.h"
#include "buffersIrregFile.h"
#include "buffersRegFile.h"
#include "json.h"
namespace SEP {
/*!
  Class for file broken into many pieces
  */
class fileBuffersRegFile : public buffersRegFile {
 public:
  /*!
    Create default object
    */
  fileBuffersRegFile() { ; }
  /*!
    Create file object

    \param arg JSON argument list
    \param usage Usage for file
    \param tag Tag for file
    \param progName Name of probram
    \param ndimMax Minimum number of dimensions for the dataset
    */
  fileBuffersRegFile(const Json::Value &arg, const SEP::usage_code usage,
                     const std::string &tag, const std::string &progName,
                     const int ndimMax = -1);
  /*!
  Remove description and directory
  */
  virtual void removeDescDir() override;
  /*! Close the file */
  virtual void close() override;
  /*! Create buffers */
  void createBuffers() override;
};
/*!
  Create file buffer for irregular file
  */
class fileBuffersIrregFile : public buffersIrregFile {
 public:
  /*!
     Create default file oject
     */
  fileBuffersIrregFile() { ; }
  /*!
  Create file object

  \param arg JSON argument list
  \param usage Usage for file
  \param tag Tag for file
  \param progName Name of probram
  \param ndimMax Minimum number of dimensions for the dataset
  */
  fileBuffersIrregFile(const Json::Value &arg, const SEP::usage_code usage,
                       const std::string &tag, const std::string &progName,
                       const int ndimMax = -1);
  /*! Close the file */

  virtual void close();
  /*! Create buffers */

  void createBuffers();

};  // namespace SEP

}  // namespace SEP

#endif
