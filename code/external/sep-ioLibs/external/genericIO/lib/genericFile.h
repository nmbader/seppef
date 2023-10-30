#ifndef GENERIC_FILE
#define GENERIC_FILE 1
#include "byteHyper.h"
#include "complexDoubleHyper.h"
#include "complexHyper.h"
#include <complex>
#include <memory>
#include <vector>

#include "doubleHyper.h"
#include "floatHyper.h"
#include "header.h"
#include "hypercube.h"
#include "intHyper.h"
#include "ioConfig.h"
#include "ioTypes.h"
#include "json.h"
#include "paramObj.h"
#include "sepVectorConfig.h"
namespace SEP {
enum usage_code {
  usageIn,    ///< File only for read operations
  usageOut,   ///< File only for write operations
  usageInOut, ///< File for both input and output
  usageScr    ///< Scratch file (remove at end)
};
/*!
   Generic regular file object (virtual)
*/
class genericRegFile : public paramObj {
public:
  /*!
   Default initialization
*/
  genericRegFile() { _type = DATA_UNKNOWN; }
  /*!
     Put integer to file

     \param par Name of parameter
     \param val Value to write to file description
  */
  /*!
    Remove file from system
    */
  virtual void remove() { perror("must override remove"); }
  virtual void putInt(const std::string &par, const int val) = 0;
  /*!
   Put float to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putFloat(const std::string &par, const float val) = 0;
  /*!
   Put string to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putString(const std::string &par, const std::string &val) = 0;
  /*!
   Put boolean to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putBool(const std::string &par, const bool val) = 0;
  /*!
   Put ints to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putInts(const std::string &par, const std::vector<int> &val) = 0;
  /*!
   Put floats to file description

   \param par Name of parameter
   \param val Value to write to file description
*/
  virtual void putFloats(const std::string &par,
                         const std::vector<float> &val) = 0;

  /*!
      Completely read a file (any type)
  */
  std::shared_ptr<regSpace> read();
  /*!
   Write entire file

   \param hyp FloatHyper (from sepVector) to store file contents into
*/
  bool writeFloatStream(std::shared_ptr<SEP::floatHyper> hyp);
  /*!
   Read entire file

   \param hyp FloatHyper (from sepVector) to store file contents into
*/
  bool readFloatStream(std::shared_ptr<SEP::floatHyper> hyp);

  /*!
 Read a portion of file based on window parameters

 \param nw,fw,jw Standard window parameters
 \param hyp Floathyper (from sepVector) storage
*/
  bool readFloatWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                       const std::vector<int> &jw,
                       std::shared_ptr<SEP::floatHyper> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp Floathyper (from sepVector) storage
*/
  bool writeFloatWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                        const std::vector<int> &jw,
                        std::shared_ptr<SEP::floatHyper> hyp);
  /*!
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  bool readByteStream(std::shared_ptr<SEP::byteHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  bool writeByteStream(const std::shared_ptr<SEP::byteHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  bool readByteWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                      const std::vector<int> &jw,
                      std::shared_ptr<SEP::byteHyper> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  bool writeByteWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                       const std::vector<int> &jw,
                       std::shared_ptr<SEP::byteHyper> hyp);

  /*!
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  bool readIntStream(std::shared_ptr<SEP::intHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  bool writeIntStream(const std::shared_ptr<SEP::intHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  bool readIntWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                     const std::vector<int> &jw,
                     const std::shared_ptr<SEP::intHyper> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp intHyper (from sepVector) storage
*/
  bool writeIntWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                      const std::vector<int> &jw,
                      std::shared_ptr<SEP::intHyper> hyp);

  /*!
Read entire file

\param hyp complexHyper (from sepVector) to grab file contents from
*/
  bool readComplexStream(std::shared_ptr<SEP::complexHyper> hyp);

  /*!
Read entire file

\param hyp complexDoubleHyper (from sepVector) to grab file contents from
*/
  bool readComplexDoubleStream(std::shared_ptr<SEP::complexDoubleHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  bool writeComplexStream(const std::shared_ptr<SEP::complexHyper> hyp);

  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  bool
  writeComplexDoubleStream(const std::shared_ptr<SEP::complexDoubleHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  bool readComplexWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                         const std::vector<int> &jw,
                         std::shared_ptr<SEP::complexHyper> hyp);

  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  bool readComplexDoubleWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               std::shared_ptr<SEP::complexDoubleHyper> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  bool writeComplexWindow(const std::vector<int> &nw,
                          const std::vector<int> &fw,
                          const std::vector<int> &jw,
                          std::shared_ptr<SEP::complexHyper> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  bool writeComplexDoubleWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                std::shared_ptr<SEP::complexDoubleHyper> hyp);

  /*!
Read entire file

\param hyp doubleHyper (from sepVector) to grab file contents from
*/
  bool readDoubleStream(std::shared_ptr<SEP::doubleHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  bool writeDoubleStream(const std::shared_ptr<SEP::doubleHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  bool readDoubleWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                        const std::vector<int> &jw,
                        const std::shared_ptr<SEP::doubleHyper> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  bool writeDoubleWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                         const std::vector<int> &jw,
                         std::shared_ptr<SEP::doubleHyper> hyp);
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp generic sepVector storage
*/
  bool readWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                  const std::vector<int> &jw,
                  std::shared_ptr<SEP::regSpace> hyp);
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  bool writeWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                   const std::vector<int> &jw,
                   std::shared_ptr<SEP::regSpace> hyp);
  /*!
Read a byte stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void readByteStream(unsigned char *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readByteStream is undefined"));
  }
  /*!
Read a float stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void readFloatStream(float *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readFloatStream is undefined"));
  }
  /*!
Write a float stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void writeFloatStream(const float *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeFloatStream is undefined"));
  }
  /*!
Write a byte stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void writeByteStream(const unsigned char *array,
                               const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeByteStream is undefined"));
  }
  /*!
Read a byte window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readByteWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw,
                              unsigned char *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("readByteWindow is undefined"));
  }
  /*! Seek to a given position inot a file
  \param iv Relative location
  \param whence (0, begining; 1, current; 2, end )
  */
  virtual void seekTo(const long long iv, const int whence) {
    if (whence == iv) {
      ;
    }
  }
  /*!
Read a complex window

\param npts Number of elements
\param array Array to read into
*/
  virtual void readComplexStream(std::complex<float> *array,
                                 const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readComplexStream is undefined"));
  }
  /*!
Read a complexDouble window

\param npts Number of elements
\param array Array to read into
*/
  virtual void readComplexDoubleStream(std::complex<double> *array,
                                       const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readComplexDoubleStream is undefined"));
  }

  /*!
Write a complex window

\param npts Number of elements
\param array Array to read into
*/
  virtual void writeComplexStream(const std::complex<float> *array,
                                  const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeComplexStream is undefined"));
  }
  /*!
Write a complex window

\param npts Number of elements
\param array Array to read into
*/
  virtual void writeComplexDoubleStream(const std::complex<double> *array,
                                        const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeComplexDoubleStream is undefined"));
  }
  /*!
Read a byte window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readComplexWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<float> *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("readComplexWindow is undefined"));
  }

  /*!
Read a byte window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readComplexDoubleWindow(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw,
                                       std::complex<double> *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("readComplexDoubleWindow is undefined"));
  }
  /*!
Write a complex window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeComplexWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<float> *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("writeComplexWindow is undefined"));
  }

  /*!
Write a complexDouble window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeComplexDoubleWindow(const std::vector<int> &nw,
                                        const std::vector<int> &fw,
                                        const std::vector<int> &jw,
                                        const std::complex<double> *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("writeComplexDoubleWindow is undefined"));
  }
  /*!
     Get the size of a dataset
  */

  virtual long long getDataSize() {
    throw SEPException(std::string("getDataSize is undefined"));
  }
  /*!
Read a float window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readFloatWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw, float *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("readFloatWindow is undefined"));
  }
  /*!
Write a float window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeFloatWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                const float *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("writeFloatWindow is undefined"));
  }
  /*!
Read a double window

\param npts Number of points
\param array Array to read into
*/
  virtual void readDoubleStream(double *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readDoubleStream is undefined"));
  }
  /*!
Write a double stream

\param npts Number of points
\param array Array to read into
*/
  virtual void writeDoubleStream(const double *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeDoubleStream is undefined"));
  }
  /*!
Read a double window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readDoubleWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw, double *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("readDoubleWindow is undefined"));
  }
  /*!
Write a double window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 const double *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("writeDoubleWindow is undefined"));
  }
  /*!
Read a integer stream

\param npts Number of points
\param array Array to read into
*/
  virtual void readIntStream(int *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readIntStream is undefined"));
  }
  /*!
Write a integer stream

\param npts Number of points
\param array Array to read into
*/
  virtual void writeIntStream(const int *array, const long long npts) {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeIntStream is undefined"));
  }
  /*!
Read a integer window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void readIntWindow(const std::vector<int> &nw,
                             const std::vector<int> &fw,
                             const std::vector<int> &jw, int *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("readIntWindow is undefined"));
  }
  /*!
Write a integer window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeIntWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw, const int *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("writeIntWindow is undefined"));
  }
  /*!
Write a byte window

\param nw,fw, jw Window parameters
\param array Array to read into
*/
  virtual void writeByteWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               const unsigned char *array) {
    if (nw.size() == 0 && fw.size() == 0 && jw.size() == 0 && array != 0)
      ;
    throw SEPException(std::string("writeByteWindow is undefined"));
  }
  /*!
    Read the description of the file
*/
  virtual void readDescription(const int ndim = -1) = 0;

  /*!
Write  the description of the file
*/
  virtual void writeDescription() { ; }
  /*!
  Close a file
  */
  virtual void close() { ; }
  /*!
    Set the hypercube for the file
    \param hyp Hypercube describing regular file
    */
  void setHyper(const std::shared_ptr<SEP::hypercube> hyp);
  /*!
    Return datatype of file
    */
  dataType getDataType() { return _type; }
  /*!
   Get the size of each element (float-4,complex-8,etc)
  */
  int getDataEsize();
  /*!
    Set the data type

    \param type Datatype
  */
  void setDataType(const dataType typ) { _type = typ; }
  /*!
    Set the data type
     \param typ (String describing datatype)
     */
  void setDataType(const std::string &typ) {
    setDataType(SEP::toElementType(typ));
  }
  /*!
    Get the data type as a string
    */
  std::string getDataTypeString();

  /*!
     Get usage for file

  */
  usage_code getUsageCode() { return _usage; }

  virtual Json::Value getDescription() = 0;

  virtual void putDescription(const std::string &title,
                              const Json::Value &desc) = 0;

  /*!
   Return hypercube describing dataset
   */
  const std::shared_ptr<SEP::hypercube> getHyper() {
    if (_hyper == nullptr)
      throw SEPException(std::string("hypercube not defined"));

    return _hyper;
  }

  virtual ~genericRegFile() { ; }

  /*!
  Return binary location

  @return binary location
  */

  std::string getBinary() const { return _binary; }

protected:
  std::string _binary;

  std::shared_ptr<SEP::hypercube> _hyper =
      nullptr;                        ///< Hypercube describing the RSF
  dataType _type = SEP::DATA_UNKNOWN; ///< The dataype for for the RSF
  usage_code _usage;
};

class genericHeaderObj {
public:
  genericHeaderObj() { ; }

  virtual void readFloatData(std::shared_ptr<SEP::floatHyper> buf) = 0;
  virtual void readFloatData(std::shared_ptr<SEP::floatHyper> buf,
                             const std::vector<bool> ind) = 0;
  virtual void writeFloatData(const std::shared_ptr<SEP::floatHyper> buf) = 0;

  virtual void readByteData(std::shared_ptr<SEP::byteHyper> buf) = 0;

  virtual void readByteData(std::shared_ptr<SEP::byteHyper> buf,
                            const std::vector<bool> ind) = 0;
  virtual void writeByteData(const std::shared_ptr<SEP::byteHyper> buf) = 0;

  virtual void readDoubleData(std::shared_ptr<SEP::doubleHyper> buf) = 0;
  virtual void readDoubleData(std::shared_ptr<SEP::doubleHyper> *buf,
                              const std::vector<bool> ind) = 0;
  virtual void writeDoubleData(const std::shared_ptr<SEP::doubleHyper> buf) = 0;

  virtual void readComplexData(std::shared_ptr<SEP::complexHyper> buf) = 0;
  virtual void readComplexData(std::shared_ptr<SEP::complexHyper> buf,
                               const std::vector<bool> ind) = 0;
  virtual void
  writeComplexData(const std::shared_ptr<SEP::complexHyper> buf) = 0;
  virtual void
  readComplexDoubleData(std::shared_ptr<SEP::complexDoubleHyper> buf) = 0;
  virtual void
  readComplexDoubleData(std::shared_ptr<SEP::complexDoubleHyper> buf,
                        const std::vector<bool> ind) = 0;
  virtual void writeComplexDoubleData(
      const std::shared_ptr<SEP::complexDoubleHyper> buf) = 0;
  virtual void readIntData(std::shared_ptr<SEP::intHyper> buf) = 0;
  virtual void readIntData(std::shared_ptr<SEP::intHyper> buf,
                           const std::vector<bool> ind) = 0;
  virtual void writeIntData(const std::shared_ptr<SEP::intHyper> buf) = 0;
  std::shared_ptr<header> getHeader();
  std::shared_ptr<header> cloneHeader() { return _header->clone(); }

private:
  std::shared_ptr<header> _header;
};

class genericIrregFile : public genericRegFile {
public:
  genericIrregFile() {}

  virtual std::shared_ptr<genericHeaderObj>
  readHeaderWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                   const std::vector<int> &jw) = 0;
  virtual void writeHeaderWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::shared_ptr<genericHeaderObj> header,
                                 std::vector<bool> &exists) = 0;
};
} // namespace SEP

#endif
