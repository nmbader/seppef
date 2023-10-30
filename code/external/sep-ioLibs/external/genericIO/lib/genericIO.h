#ifndef GENERIC_IO_H
#define GENERIC_IO_H 1
#include <map>
#include <memory>
#include "genericFile.h"
#include "json.h"
#include "paramObj.h"

namespace SEP {
/*! Abstract Class for different types IO*/
class genericIO {
 public:
  /*!
 Default object, does nothing
*/
  genericIO() { _type = "NONE"; }

  /*!
     Return a genericRegFile object

    \param tag Tag used to access dataset
    \param name Name of dataset
    \param usage Usage for file (in,out,scratch)
    \param ndimMax Output file should have ndimMax axes
  */
  std::shared_ptr<SEP::genericRegFile> getRegFile(const std::string &name,
                                                  const SEP::usage_code usage,
                                                  const int ndimMax = -1);

  /*!
     Return a genericRegFile object

    \param tag Tag used to access dataset
    \param doc Documentataion for file
    \param name Name of dataset
    \param usage Usage for file (in,out,scratch)
    \param ndimMax Output file should have ndimMax axes
  */
  std::shared_ptr<SEP::genericRegFile> getDocRegFile(
      const std::string &name, const std::string &doc,
      const SEP::usage_code usage, const int ndimMax = -1);

  /*!
   Return a genericIrregFile object

  \param tag Tag used to access dataset
  \param name Name of dataset
  \param usage Usage for file (in,out,scratch)
  \param ndimMax Output file should have ndimMax axes

*/
  std::shared_ptr<SEP::genericIrregFile> getIrregFile(
      const std::string &name, const SEP::usage_code usage,
      const int ndimMax = -1);

  /*!
     Return a genericRegFile object

    \param tag Tag used to access dataset
    \param name Name of dataset
    \param usage Usage for file (in,out,scratch)
    \param ndimMax Output file should have ndimMax axes
  */
  virtual std::shared_ptr<SEP::genericRegFile> getRegFileTag(
      const std::string &tag, const std::string &name,
      const SEP::usage_code usage, const int ndimMax = -1) = 0;

  /*!
     Return a genericRegFile object

    \param name Name of dataset
    \param doc Documentation for file
    \param usage Usage for file (in,out,scratch)
    \param ndimMax Output file should have ndimMax axes
  */
  virtual std::shared_ptr<SEP::genericRegFile> getDocRegFile(
      const std::string &name, const std::string &doc, const std::string usage,
      const int ndimMax = -1);

  /*!
     Return a genericRegFile object

    \param name Name of dataset
    \param usage Usage for file (in,out,scratch)
    \param ndimMax Output file should have ndimMax axes
  */
  virtual std::shared_ptr<SEP::genericRegFile> getRegFile(
      const std::string &name, const std::string usage,
      const int ndimMax = -1) {
    SEP::usage_code code;
    if (usage == std::string("UsageIn")) {
      code = usageIn;
    } else if (usage == std::string("UsageInOut"))
      code = usageInOut;
    else if (usage == std::string("UsageOut"))
      code = usageOut;
    else if (usage == std::string("UsageScr"))
      code = usageScr;
    else {
      _param->error(std::string("Unknown code ") + usage);
    }
    return getRegFile(name, code);
  }
  /*!
   Return a genericIrregFile object

  \param tag Tag used to access dataset
  \param name Name of dataset
  \param usage Usage for file (in,out,scratch)
  \param ndimMax Output file should have ndimMax axes

*/
  virtual std::shared_ptr<SEP::genericIrregFile> getIrregFileTag(
      const std::string &tag, const std::string &name,
      const SEP::usage_code usage, const int ndimMax = -1) = 0;

  /*!
    Return parameter object associated with this IO
    */
  virtual std::shared_ptr<paramObj> getParamObj() { return _param; }
  /*!
     Add file to the list of regular files being used by this IO type

     \param fle File to add the list of files
  */
  void addRegFile(std::string name, std::shared_ptr<genericRegFile> fle) {
    _regFiles[name] = fle;
  }

  /*!
     Replace parameter type

     \param obj Replace parameter object
  */
  void replaceParamObj(std::shared_ptr<paramObj> obj) { _param = obj; }

  /*!
   Add file to the list of irregular files being used by this IO type

   \param fle File to add the list of files
*/
  void addIrregFile(std::string name, std::shared_ptr<genericIrregFile> fle) {
    _irregFiles[name] = fle;
  }
  /*!
     Return the file object associated with a given name

          \param name Tag for the file
  */

  std::shared_ptr<genericRegFile> getRegFile(const std::string name) {
    if (!regFileExists(name))
      _param->error(std::string("Requested unknown file ") + name);
    return _regFiles[name];
  }
  /*!
     Check to see if a tag exist for the given IO regular file list

     \param name Name of the tag
  */
  bool regFileExists(const std::string name) {
    if (_regFiles.count(name) == 0) return false;
    return true;
  }
  /*!
      Check to see if a tag exist for the given IO irregular file list

      \param name Name of the tag
   */
  bool irregFileExists(const std::string name) {
    if (_irregFiles.count(name) == 0) return false;
    return true;
  }
  /*!
    Write the contents of data to a file for debugging purposes

    \param name Name of the file to write
    \param data Contents to write the file
    \param n1,n2 Dimensions of the dataset
 */
  void fileDebug(const std::string name, const float *data, const int n1,
                 const int n2);
  /*!
  Write the contents of data to a file for debugging purposes

  \param name Name of the file to write
  \param data Contents to write the file
  \param n1,n2,n3 Dimensions of the dataset
*/
  void fileDebug(const std::string name, const float *data, const int n1,
                 const int n2, const int n3);
  /*!
     Return the file object associated with a given name

          \param name Tag for the file
  */
  std::shared_ptr<genericIrregFile> getIrregFile(const std::string x) {
    if (_irregFiles.count(x) == 0)
      _param->error(std::string("Requested unknown file ") + x);
    return _irregFiles[x];
  }
  /*!
     Whether or not this IO is a valid one to use

  */
  void setValid(const bool x) { _valid = x; }
  /*!
     Set whether or not IO type is valid to use
  */
  bool getValid() { return _valid; }

  /*!
   Close all files

*/
  virtual void close() { filesClose(); }

  /*!
     Copy info from one file to another

     \param fileIn Input file to grab history from

     \param title Title in ouput description to describe info

     \param fileOut  Output description to write to
     */

  void addFileDescription(const std::shared_ptr<genericRegFile> fileIn,
                          const std::string &title,
                          std::shared_ptr<genericRegFile> fileOut);

  /*!
    Remove file from the system
    */
  void removeRegFile(const std::string fle);
  /*!
 Delete IO type (close files)

*/
  ~genericIO() { close(); }

  void setParamObj(std::shared_ptr<paramObj> par) { _param = par; }
  /*!
  Close all files

*/

  /*!
     @return ioType
   */
  std::string getType() const { return _type;}
  virtual void filesClose();
  std::string _type;  ///< IO type

 protected:
  std::map<std::string, std::shared_ptr<genericRegFile> >
      _regFiles;  ///< All active tags for regular files
  std::map<std::string, std::shared_ptr<genericIrregFile> >
      _irregFiles;                   ///< All active tags for irregular files
  std::shared_ptr<paramObj> _param;  ///< IO handler for this IO type
  bool _valid;                       ///< Whether or not IO type is valid
};                                   // namespace SEP
}  // namespace SEP
#endif
