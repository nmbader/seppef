#ifndef BUFFERS_IO
#define BUFFERS_IO_H 1
#include "buffersRegFile.h"
#include "jsonGenericIO.h"
namespace SEP {
typedef int MyCustomType;
/*!
Handled buffer (multi block) UO
*/
class buffersIO : public jsonGenericIO {
 public:
  /*!
   Default object, does nothing
*/
  buffersIO() { ; }
  /*!
  Initialize IO with command line arguments

  \param argc C number of command line arguments
  \param argv List of command line arguments
*/
  buffersIO(const int argc, char **argv) { initJsonPars(argc, argv); }
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
  virtual void close();
  /*!
   Get parameter object

*/
  virtual std::shared_ptr<paramObj> getParamObj();

 private:
  bool _init = false;  ///< Whether or not IO has been initialized
};
}  // namespace SEP
#endif
