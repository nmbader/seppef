#ifndef buffers_reg_file_h
#define buffers_reg_file_h 1
#include <stdbool.h>
#include <string>
#include "basicIO.h"
#include "buffers.h"
#include "json.h"
#include "jsonGenericRegFile.h"
namespace SEP {
/*!
Object for handling regular files broken into blocks
*/
class buffersRegFile : public jsonGenericRegFile {
 public:
  /*!
     Initialize an irregular file
 */
  buffersRegFile() { ; }
  /*!
     Write the description for the file
 */
  virtual void remove() override {
    _bufs->remove();
    removeDescDir();
  }
  /*!
  Remove description and directory
  */
  virtual void removeDescDir() {
    perror("must override description and directory");
  }
  /*!
  Write description
  */
  virtual void writeDescription() override;
  /*!
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  virtual void readByteStream(unsigned char *array,
                              const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexStream(const std::complex<float> *array,
                                  const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }

  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  virtual void writeComplexDoubleStream(const std::complex<double> *array,
                                  const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }

  /*!
Read entire file

\param hyp complexHyper (from sepVector) to grab file contents from
*/
  virtual void readComplexStream(std::complex<float> *array,
                                 const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }

    /*!
Read entire file

\param hyp complexDoubleHyper (from sepVector) to grab file contents from
*/
  virtual void readComplexDoubleStream(std::complex<double> *array,
                                 const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Write a float stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void writeFloatStream(const float *array,
                                const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Read a float stream

\param array Array to read into
\param npts Number of values to read
*/
  virtual void readFloatStream(float *array, const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }

  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void writeIntStream(const int *array, const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Read entire file

\param hyp byteHyper (from sepVector) to grab file contents from
*/
  virtual void readIntStream(int *array, const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  virtual void writeDoubleStream(const double *array,
                                 const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Read entire file

\param hyp doubleHyper (from sepVector) to grab file contents from
*/
  virtual void readDoubleStream(double *array, const long long npts) override {
    error(std::string("can not stream buffer datasets, must use window"));
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void readByteWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw,
                              unsigned char *array) override {
    setDataType(DATA_FLOAT);
    createBuffers();
    _bufs->getWindow(nw, fw, jw, (void *)array);
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp Floathyper (from sepVector) storage
*/
  virtual void readFloatWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               float *array) override {
    setDataType(DATA_FLOAT);
    createBuffers();

    _bufs->getWindow(nw, fw, jw, (void *)array);
  }
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp Floathyper (from sepVector) storage
*/
  virtual void writeFloatWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                const float *array) override {
    setDataType(DATA_FLOAT);
    createBuffers();

    _bufs->putWindow(nw, fw, jw, (void *)array);
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void readComplexWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<float> *array) override {
    createBuffers();
    _bufs->getWindow(nw, fw, jw, (void *)array);
  }

    /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  virtual void readComplexDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<double> *array) override {
    createBuffers();
    _bufs->getWindow(nw, fw, jw, (void *)array);
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexHyper (from sepVector) storage
*/
  virtual void writeComplexWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<float> *array) override {
    createBuffers();
    _bufs->putWindow(nw, fw, jw, (void *)array);
  }

    /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp complexDoubleHyper (from sepVector) storage
*/
  virtual void writeComplexDoubleWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  const std::complex<double> *array) override {
    createBuffers();
    _bufs->putWindow(nw, fw, jw, (void *)array);
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp byteHyper (from sepVector) storage
*/
  virtual void readIntWindow(const std::vector<int> &nw,
                             const std::vector<int> &fw,
                             const std::vector<int> &jw, int *array) override {
    createBuffers();
    _bufs->getWindow(nw, fw, jw, (void *)array);
  }
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp intHyper (from sepVector) storage
*/
  virtual void writeIntWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw,
                              const int *array) override {
    createBuffers();
    _bufs->putWindow(nw, fw, jw, (void *)array);
  }
  /*!
Read a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  virtual void readDoubleWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                double *array) override {
    createBuffers();
    _bufs->getWindow(nw, fw, jw, (void *)array);
  }
  /*!
Write a portion of file based on window parameters

\param nw,fw,jw Standard window parameters
\param hyp doubleHyper (from sepVector) storage
*/
  virtual void writeDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 const double *array) override {
    createBuffers();
    _bufs->putWindow(nw, fw, jw, (void *)array);
  }
  /*!
    Set memory usage object
    \param mem Memory object
    */
  void setMemoryUsage(std::shared_ptr<SEP::IO::memoryUsage> mem) {
    if (!_hyper) throw SEPException(std::string("Hypercube has not been set"));
    _mem = mem;
  }
  /*!
  Set compression usage object
  \param com Compression object
  */
  void setCompression(std::shared_ptr<SEP::IO::compress> com) { _comp = com; }
  /*!
  Set block usage object
  \param block Block object
  */
  void setBlocking(std::shared_ptr<SEP::IO::blocking> block) { _block = block; }
  /*!
 Create buffers
  */
  virtual void createBuffers() = 0;

 protected:
  std::shared_ptr<SEP::IO::buffers> _bufs = nullptr;     ///< Buffer object
  std::shared_ptr<SEP::IO::memoryUsage> _mem = nullptr;  ///< Memory object
  std::shared_ptr<SEP::IO::compress> _comp = nullptr;    ///< Compression object
  std::shared_ptr<SEP::IO::blocking> _block = nullptr;   ///< Blocking object

};  // namespace SEP

}  // namespace SEP
#endif
