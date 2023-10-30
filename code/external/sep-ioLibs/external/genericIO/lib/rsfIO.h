#ifndef RSF_IO_H
#define RSF_IO_H 1
#include "genericIO.h"
#include "rsfParams.h"
#include "rsfRegFiles.h"
/*!
  Class for handling IO operations when writing to multiple files/objects
*/
class rsfIO : public genericIO {
 public:
  /*!
Initialize IO with command line arguments

\param argc C number of command line arguments
\param argv List of command line arguments
*/
  rsfIO(int argc, char **argv);
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
      const int ndimMax = -1) override {
    throw SEPException(std::string("Undefined getIrregFileTag"));
  }
};
#endif