#ifndef SEGY_H
#define SEGY_H 1
#include "jsonGenericIO.h"
namespace SEP {
class segyIO : public jsonGenericIO {
 public:
  /*!
Initialize IO with command line arguments

\param argc C number of command line arguments
\param argv List of command line arguments
*/
  segyIO(const int argc, char **argv) { initJsonPars(argc, argv); 
  _type="SEGY";
}
  /*!
   Return a genericRegFile object

  \param tag Tag used to access dataset
  \param name Name of dataset
  \param usage Usage for file (in,out,scratch)
  \param ndimMax Output file should have ndimMax axes
*/
  virtual std::shared_ptr<genericRegFile> getRegFileTag(const std::string &tag,
                                                        const std::string &name,
                                                        const usage_code usage);
  /*!
 Return a genericIrregFile object

\param tag Tag used to access dataset
\param name Name of dataset
\param usage Usage for file (in,out,scratch)
\param ndimMax Output file should have ndimMax axes

*/
  virtual std::shared_ptr<genericIrregFile> getIrregFileTag(
      const std::string &tag, const std::string &name, const usage_code usage) {
    throw(SEPException(std::string("getIrregFileTag undefined for segyIO")));
  }
  /*!
Delete SEGYIO object

*/
  ~segyIO() { ; }
};

}  // namespace SEP
#endif
