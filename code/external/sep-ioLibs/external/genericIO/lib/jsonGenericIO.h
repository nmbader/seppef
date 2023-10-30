#ifndef JSON_IO_H
#define JSON_IO_H 1
#include <fstream>  // std::ifstream
#include "genericIO.h"
#include "jsonGenericIrregFile.h"
#include "jsonGenericRegFile.h"

#include "jsonParamObj.h"
namespace SEP {
/*!
Class for handling IO for files written using JSON
*/

class jsonGenericIO : public genericIO {
 public:
  /*!
 Default object, does nothing
*/
  jsonGenericIO() { ; }
  /*!
Initialize IO with command line arguments

\param argc C number of command line arguments
\param argv List of command line arguments
*/
  jsonGenericIO(const int argc, char **argv) { initJsonPars(argc, argv); }
  /*!
Initialize IO JSON parameters
\param argc C number of command line arguments
\param argv List of command line arguments
*/
  void initJsonPars(const int argc, char **argv);
  /*!
     Return JSON arguments
  */
  Json::Value getArgs() { return jsonArgs; }
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
 Get parameter object

*/
  virtual std::shared_ptr<paramObj> getParamObj() override;

 protected:
  std::string _progName;  ///< Name of the program
  Json::Value jsonArgs;   ///< JSON structure containing arguments
  bool _sentError =
      false;  ///< Whether or not an error on initialization has been sent

 private:
  std::ifstream inps;  ///< Stream object for IO
  bool _init;          ///< Whether or not IO has been initialized
};
}  // namespace SEP
#endif
