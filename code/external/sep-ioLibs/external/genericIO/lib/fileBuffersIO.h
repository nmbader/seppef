#ifndef file_buffers_io_h
#define file_buffers_io_h 1
#include "fileBuffersRegFile.h"
#include "jsonGenericIO.h"
namespace SEP {
/*!
  Class for handling buffers (multiple blocks) writen to conventional disk
*/
class fileBuffersIO : public jsonGenericIO {
 public:
  /*!
   Default object, does nothing
*/
  fileBuffersIO() { ; }
  /*!
Initialize IO with command line arguments

\param argc C number of command line arguments
\param argv List of command line arguments
*/
  fileBuffersIO(const int argc, char **argv) { initJsonPars(argc, argv); }
  /*!
     Return a genericRegFile object

    \param tag Tag used to access dataset
    \param name Name of dataset
    \param usage Usage for file (in,out,scratch)
    \param ndimMax Output file should have ndimMax axes
  */
  virtual std::shared_ptr<genericRegFile> getRegFileTag(
      const std::string &tag, const std::string &name, const usage_code usage,
      const int ndimMax = -1) override;
  /*!
   Return a genericIrregFile object

  \param tag Tag used to access dataset
  \param name Name of dataset
  \param usage Usage for file (in,out,scratch)
  \param ndimMax Output file should have ndimMax axes

*/

  virtual std::shared_ptr<genericIrregFile> getIrregFileTag(
      const std::string &tag, const std::string &name, const usage_code usage,
      const int ndimMax = -1) override;
  /*!
    Close all files

*/
  virtual void close() override;
  /*!
  Return parameter object associated with this IO
  */
  virtual std::shared_ptr<paramObj> getParamObj() override;

 private:
  bool _init = false;  ///< Whether or not IO has been initailized
};
}  // namespace SEP
#endif
