#ifndef JSON_IRREGFILE_FUNC_H
#define JSON_IRREGFILE_FUNC_H 1
#include <stdbool.h>
#include <string>
#include "basicIO.h"
#include "genericFile.h"
#include "ioTypes.h"
#include "json.h"
namespace SEP {

/*!
   Irregular file described by JSON
*/
class jsonGenericIrregFile : public genericIrregFile {
 public:
  // sepRegFile::sepRegFile(const std::string tag,usage_code usage){
  /*!
   Default initialization
*/
  jsonGenericIrregFile() { ; }
  /*!
    \param arg JSON arguments
    \param usage Usage for file
    \param traceH Length of trace headers
    \param progName Name of program
    \param Minimum number of dimensions for dataset
    */
  jsonGenericIrregFile(const Json::Value &arg, const SEP::usage_code usage,
                       const std::string &tag, const int reelH,
                       const int traceH, const std::string &progName,
                       const int ndim = -1);
  /*!
     Setup JSON environmnet
    \param jsonArgs Arguments described by JSON
    \param tag Tag for dataset
    \param desFileDefault Default name for output description file
  */
  void setupJson(const Json::Value &jsonArgs, const std::string &tag,
                 const std::string desFileDefault = std::string(""));
  /*!
   Get an integer from a file

   \param arg Name of the prameter
  */
  virtual int getInt(const std::string &arg) const override;
  /*!
  Get an integer from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
 */
  virtual int getInt(const std::string &arg, const int def) const override;
  /*!
  Get a float from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
 */
  virtual float getFloat(const std::string &, const float def) const override;
  /*!
Get a float from a file

\param arg Name of the prameter
*/
  virtual float getFloat(const std::string &) const override;
  /*!
  Get a string  from a file

  \param arg Name of the prameter
 */
  virtual std::string getString(const std::string &arg) const override;
  /*!
Get a string from a file

\param tag Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::string getString(const std::string &arg,
                                const std::string &def) const override;
  /*!
  Get boolean from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
 */
  virtual bool getBool(const std::string &, const bool def) const override;
  /*!
Get a boolean from a file

\param arg Name of the prameter
*/
  virtual bool getBool(const std::string &) const override;
  /*!
  Get integer from a file

  \param arg Name of the prameter
  \param nvals Number of values
 */
  virtual std::vector<int> getInts(const std::string &arg,
                                   int nvals) const override;
  /*!
Get integers from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::vector<int> getInts(const std::string &arg,
                                   const std::vector<int> &defs) const override;
  /*!
  Get an floats from a file

  \param arg Name of the prameter
  \param nval Number of values to look for
 */
  virtual std::vector<float> getFloats(const std::string &arg,
                                       const int nvals) const override;
  /*!
Get floats from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::vector<float> getFloats(
      const std::string &arg, const std::vector<float> &defs) const override;
  /*!
    Output a message
    \param err Message to output
    */
  virtual void message(const std::string &err) const override;
  /*!
    Output a message and exit with an error
    \param err Message to output
    */
  virtual void error(const std::string &err) const override;
  /*!
      Return JSON arguments
  */
  Json::Value getArgs() { return jsonArgs; }
  /*!
    Get tag associated with file
    */
  std::string getTag() { return _tag; }
  /*!
  Close a file
  */
  virtual void close() override;
  /*!
    Set usage of file
    \param usage Usage for file
  */
  void setUsage(const usage_code usage) { _usage = usage; }
  /*!
     Set the history file file
     */
  void setHistory(const Json::Value &hist);
  /*!
     Get usage for file
     */
  usage_code getUsage() { return _usage; }
  /*!
  Get the name of the JSON file
  */
  virtual std::string getJSONFileName() const;
  /*!
      Get the name of the data file
      */
  std::string getDataFileName() const;
  /*!
     Read file description, with minimum nuber of dimensions
     */
  virtual void readDescription(const int ndimMax=-1) override;
  /*!
     Write file description
     */
  virtual void writeDescription() override;
  /*!
   Put integer to file

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putInt(const std::string &par, const int val) override;
  /*!
   Put float to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putFloat(const std::string &par, const float val) override;
  /*!
 Put string to file description

 \param par Name of parameter
 \param val Value to write to file description
*/
  virtual void putString(const std::string &par,
                         const std::string &val) override;
  /*!
   Put boolean to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putBool(const std::string &par, const bool val) override;
  /*!
 Put ints to file description

 \param par Name of parameter
 \param val Value to write to file description
*/
  virtual void putInts(const std::string &par,
                       const std::vector<int> &val) override;
  /*!
   Put floats to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putFloats(const std::string &par,
                         const std::vector<float> &val) override;
  virtual std::shared_ptr<genericHeaderObj> readHeaderWindow(
      const std::vector<int> &nw, const std::vector<int> &fw,
      const std::vector<int> &jw) override {}
  virtual void writeHeaderWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::shared_ptr<genericHeaderObj> header,
                                 std::vector<bool> &exists) override {
    ;
  }

  virtual Json::Value getDescription() override { return jsonArgs; }
  /*!

    \param title Name to give the history

    \param desc Description to putDescription
*/

  virtual void putDescription(const std::string &title,
                              const Json::Value &desc) override {
    jsonArgs[title] = desc;
  }

 protected:
  Json::Value jsonArgs;  ///< JSON values
  bool _newFile;         ///< Whether or not a new file
  dataType _dtype;       ///< Datatype associated with the dataset

  std::string _tag;                ///< Tag assicated with the dataset
  usage_code _usage;               ///< Usage file
  std::string _jsonFile;           ///< Description file
  std::string _dataFile;           ///< Binary data file
  std::shared_ptr<myFileIO> myio;  ///< Basic IO ovject
  int _reelH;                      ///< Reel header
  int _traceH;                     ///< Trace header
};                                 // namespace SEP
}  // namespace SEP

#endif
