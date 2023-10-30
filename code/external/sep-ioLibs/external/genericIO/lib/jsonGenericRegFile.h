#ifndef JSON_REGFILE_FUNC_H
#define JSON_REGFILE_FUNC_H 1
#include <stdbool.h>
#include <string>
#include "basicIO.h"
#include "genericFile.h"
#include "ioTypes.h"
#include "json.h"
namespace SEP {
/*!
   Regular file described by JSON
*/
class jsonGenericRegFile : public genericRegFile {
 public:
  /*!
   Default initialization
*/
  jsonGenericRegFile() { ; }
  /*!
  \param arg JSON arguments
  \param usage Usage for file
  \param traceH Length of trace headers
  \param progName Name of program
  \param Minimum number of dimensions for dataset
  */
  jsonGenericRegFile(const Json::Value &arg, const SEP::usage_code usage,
                     const std::string &tag, const int reelH, const int traceH,
                     const std::string &progName, const int ndim = -1);
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
Get the name of the JSON file, no path
*/
  virtual std::string getJSONFileBaseName() const;



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
*/ virtual void putFloats(const std::string &par,
                          const std::vector<float> &val) override;

  /*!
   Read entire file

   \param hyp FloatHyper (from sepVector) to store file contents into
*/
  virtual void readFloatStream(float *array, const long long npts) override;
  /*!
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/

  virtual void readByteStream(unsigned char *array,
                              const long long npts) override;

  /*! Seek to a given position inot a file
   \param iv Relative location
   \param whence (0, begining; 1, current; 2, end )
   */
  virtual void seekTo(const long long iv, const int whence) override;

  /*!
 Write entire file

 \param hyp FloatHyper (from sepVector) to grab file contents from
*/
  virtual void writeFloatStream(const float *array,
                                const long long npts) override;
  /*!
Write entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  virtual void writeByteStream(const unsigned char *array,
                               const long long npts) override;
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void writeByteWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               unsigned char const *array) override;
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void readByteWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw,
                              unsigned char *array) override;
  /*!
Read a float window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readFloatWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               float *array) override;
  /*!
   Get the size of a dataset
*/
  virtual long long getDataSize() override;
  /*!
Write a float window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeFloatWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                const float *array) override;
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void readComplexWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<float> *array) override;
 
   /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  virtual void readComplexDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<double> *array) override;
 
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<float> *array) override;
   /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  virtual void writeComplexDoubleWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<double> *array) override;
 
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexStream(const std::complex<float> *array,
                                  const long long npts) override;
 
   /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexDoubleStream(const std::complex<double> *array,
                                  const long long npts) override;
 
  /*!
Read entire file

\param hyp complexHyper (from sepVector) to grab file contents from
*/
  virtual void readComplexStream(std::complex<float> *array,
                                 const long long npts) override;
 
  /*!
Read entire file

\param hyp complexDoubleHyper (from sepVector) to grab file contents from
*/
  virtual void readComplexDoubleStream(std::complex<double> *array,
                                 const long long npts) override;

  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  virtual void writeDoubleStream(const double *array,
                                 const long long npts) override;
  /*!
Read entire file

\param hyp doubleHyper (from sepVector) to grab file contents from
*/
  virtual void readDoubleStream(double *array, const long long npts) override;
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  virtual void readDoubleWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                double *array) override;
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  virtual void writeDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 const double *array) override;
  /*!
Write a integer stream

\param npts Number of points
\param array Array to read into
*/
  virtual void writeIntStream(const int *array, const long long npts) override;
  /*!
Read a integer stream

\param npts Number of points
\param array Array to read into
*/
  virtual void readIntStream(int *array, const long long npts) override;
  /*!
Read a integer window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readIntWindow(const std::vector<int> &nw,
                             const std::vector<int> &fw,
                             const std::vector<int> &jw, int *array) override;
  /*!
Write a integer window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeIntWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw,
                              const int *array) override;
  /*!
    Get the description of pthe current file

    */
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
};
}  // namespace SEP

#endif
