#ifndef MEMORY_REGFILE_FUNC_H
#define MEMORY_REGFILE_FUNC_H 1
#include <stdbool.h>
#include <string>
#include "genericFile.h"
namespace SEP {
/*!
  Regular file object for in memory (more of a testing class)

*/
class memoryRegFile : public SEP::genericRegFile {
 public:
  /*!
      Create a regular file object
       \param tg Tag for dataset
       \param usage Usage for file
       \param ndimMax Minimum number of dimensions for file
  */
  memoryRegFile(const std::string &tg, const SEP::usage_code usage,
                const int ndimMax = -1);
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
  /*! Seek to a given position inot a file
  \param iv Relative location
  \param whence (0, begining; 1, current; 2, end )
  */
  virtual void seekTo(const long long iv, const int whence) override;
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
                                   int nvals = 1) const override;
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
                                       int nvals = 1) const override;
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
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  virtual void readByteStream(unsigned char *array,
                              const long long npts) override;
  /*!
Write entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  virtual void writeByteStream(const unsigned char *array,
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
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexStream(const std::complex<float> *array,
                                  const long long npts) override;
  
    /*!
Read a portion of file based on window parameters

\param array Array to write out
\param npts Number of points to write
*/
  virtual void writeComplexDoubleStream(const std::complex<double> *array,
                                  const long long npts) override;
  
  /*!
Read a float stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void readFloatStream(float *array, const long long npts) override;
  /*!
Write a float stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void writeFloatStream(const float *array,
                                const long long npts) override;
  /*!
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  virtual void readIntStream(int *array, const long long npts) override;
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void writeIntStream(const int *array, const long long npts) override;
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
  virtual void writeDoubleStream(const double *array,
                                 const long long npts) override;
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
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void readIntWindow(const std::vector<int> &nw,
                             const std::vector<int> &fw,
                             const std::vector<int> &jw, int *array) override;
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp intHyper (from sepVector) storage
*/
  virtual void writeIntWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw,
                              const int *array) override;
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
 Read a portion of file based on window parameters

 \param nw,fw,jw Standard window parameters
 \param hyp Floathyper (from sepVector) storage
*/
  virtual void readFloatWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               float *array) override;
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void writeByteWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               const unsigned char *array) override;
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp Floathyper (from sepVector) storage
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
\param hyp complexHyper (from sepVector) storage
*/
  virtual void readComplexDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<double> *array) override;
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<float> *array) override;


  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexDoubleWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<double> *array) override;
  /*!
  Close a file
  */
  virtual void close() override;
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

  /*!
   \param str
   \return list split ,

*/
  std::vector<std::string> splitString(const std::string &str) const;

  /*!
    Read the description of the file
*/
  virtual void readDescription(const int ndim=-1) override { ; };
  virtual Json::Value getDescription() override { 
Json::Value a;
return a; }
  virtual void putDescription(const std::string &title,
                              const Json::Value &desc) override {
    ;
  }
  /*!
Write  the description of the file
*/
  virtual void writeDescription() override { ; }

  /*!
  Get storage pointer

  \return pointer to memory
  */
  void *getPtr() { return (void *)_buf.data(); }

  /*! get Vector
   \return vector
   */
  std::vector<unsigned char> returnVec() { return _buf; }

  /*!
  Allocate if not allocated
  */
  void allocateCheck(dataType typ);

 private:
  std::string _tag;    ///< Tag associated with dataset
  long long _pos = 0;  ///< Position in file
  std::map<std::string, std::string>
      _dict;  ///< Dictionary Containing parameters
  std::vector<unsigned char> _buf;
};
}  // namespace SEP

#endif
