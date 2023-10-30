#ifndef rsf_REGFILE_FUNC_H
#define rsf_REGFILE_FUNC_H 1
#include <stdbool.h>
#include <string>
#include "SEPException.h"
#include "basicIO.h"
#include "genericFile.h"
#define NO_BLAS 1
extern "C" {
#include <rsf.h>
}
/*!
  Object for handling regular file IO
*/
class rsfBasic : public basicIO {
 public:
  /*!
    Initialize RSF IO
    \param file File object
    */
  rsfBasic(sf_file file) { _file = file; }
  /*! Seek to a given position inot a file
   \param pos Relative location
   */
  virtual inline void seekToPos(const long long pos) {
    long long ft = ftell(_myf);
    long long bg = 1024 * 1024 * 1024;
    long long diff = pos - ft;
    while (diff != 0) {
      long long dst;
      if (diff > 0)
        dst = std::min(bg, diff);
      else
        dst = -std::min(-diff, bg);
      sf_seek(_file, dst, SEEK_CUR);
      diff -= dst;
    }
  }
  /*!
    Read a data stream
    \param sz Number of bytes
    \param data Storage
    */
  virtual void readStream(const long long sz, void *data) {
    sf_ucharread((unsigned char *)data, sz, _file);
  }
  /*!
   Write a data stream
   \param sz Number of bytes
   \param data Storage
   */
  virtual void writeStream(const long long sz, const void *data) {
    sf_ucharwrite((unsigned char *)data, sz, _file);
  }

  sf_file _file;  ///< RSF object
};
/*!
   RSF implementation of a generic file
   */

class rsfRegFile : public genericRegFile {
 public:
  /*!
    Initialize RSF reg file
     \param tag Tag of dataset
     \param usage Usage for file
     */
  rsfRegFile(std::string tg, usage_code usage);
  virtual int getInt(const std::string arg) const override;
  virtual int getInt(const std::string arg, const int def) const override;

  virtual float getFloat(const std::string, const float def) const override;
  virtual float getFloat(const std::string) const override;

  virtual std::string getString(const std::string arg) const override;
  virtual std::string getString(const std::string arg,
                                const std::string def) const override;

  virtual bool getBool(const std::string, const bool def) const override;
  virtual bool getBool(const std::string) const override;
  virtual void seekTo(const long long iv, const int whence) override;

  virtual std::vector<int> getInts(const std::string arg,
                                   int num) const override;
  virtual std::vector<int> getInts(const std::string arg,
                                   std::vector<int> &defs) const override;

  virtual std::vector<float> getFloats(const std::string arg,
                                       int num) const override;
  virtual std::vector<float> getFloats(const std::string arg,
                                       std::vector<float> &defs) const;

  virtual void error(const std::string err) const;
  virtual void readComplexStream(std::complex<float> *array,
                                 const long long npts);
  virtual void readFloatStream(float *array, const long long npts);
  virtual void readByteStream(unsigned char *array, const long long npts);
  virtual void writeComplexStream(const std::complex<float> *array,
                                  const long long npts);
  virtual void readIntStream(int *array, const long long npts) override {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readIntStream not defined for RSF"));
  }
  virtual void readDoubleStream(double *array, const long long npts) override {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("readDoubleStream not defined for RSF"));
  }
  virtual void writeFloatStream(const float *array,
                                const long long npts) override;
  virtual void writeIntStream(const float *array,
                              const long long npts) override {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeIntStream not defined for RSF"));
  }
  virtual void writeDoubleStream(const double *array,
                                 const long long npts) override {
    if (array == 0 && npts == 0)
      ;
    throw SEPException(std::string("writeDoubleStream not defined for RSF"));
  }

  virtual void readByteWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw, unsigned char *array);

  virtual void readComplexWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 std::complex<float> *array);

  virtual void writeComplexWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw,
                                  std::complex<float> *array);

  virtual void readFloatWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw,
                               float *array) override;

  virtual void readDoubleWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                double *array) override {
    if (fw.size() == 0 && jw.size() == 0 && nw.size() == 0 && array == 0)
      ;
    throw SEPException(std::string("readDoubleWindow not defined for RSF"));
  }
  virtual void readIntWindow(const std::vector<int> &nw,
                             const std::vector<int> &fw,
                             const std::vector<int> &jw, int *array) override {
    if (fw.size() == 0 && jw.size() == 0 && nw.size() == 0 && array == 0)
      ;
    throw SEPException(std::string("readIntWindow not defined for RSF"));
  }

  virtual void writeFloatWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw,
                                float *array) override;
  virtual void writeDoubleWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 double *array) override {
    if (fw.size() == 0 && jw.size() == 0 && nw.size() == 0 && array == 0)
      ;
    throw SEPException(std::string("writeDoubleWindow not defined for RSF"));
  }
  virtual void writeIntWindow(const std::vector<int> &nw,
                              const std::vector<int> &fw,
                              const std::vector<int> &jw, int *array) override {
    if (fw.size() == 0 && jw.size() == 0 && nw.size() == 0 && array == 0)
      ;
    throw SEPException(std::string("writeIntWindow not defined for RSF"));
  }
  virtual void readDescription(int ndimMax = -1);
  virtual void writeDescription();
  virtual void putInt(const std::string par, const int val);
  virtual void putFloat(const std::string par, const float val);
  virtual void putString(const std::string par, const std::string val);
  virtual void putBool(const std::string par, const bool val);
  virtual void putInts(const std::string par, const std::vector<int> val);
  virtual void putFloats(const std::string par, const std::vector<float> val);
  virtual void close() const { sf_fileclose(_file); }

 private:
  usage_code _usage;
  std::string _tag;
  sf_file _file;
  std::shared_ptr<rsfBasic> myio;
  int _esize;
};

#endif
