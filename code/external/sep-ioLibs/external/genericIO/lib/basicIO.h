#ifndef BASIC_IOH
#define BASIC_IOH 1
#include <stdio.h>
#include <memory>
#include "genericFile.h"
namespace SEP {
/*!
 this is very rudimentary class for file io */

class basicIO {
 public:
  /*!
    Create a basic IO classs
    */
  basicIO() { _swapData = false; }
  /*!
    Read a window
     \param nw,fw,jw  Windowing parmeters
     \param data Data associated with window
     \param head Header sassoicated with window
*/
  void readWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                  const std::vector<int> &jw, void *data, void *head = 0);
  /*!
  Write a window
   \param nw,fw,jw  Windowing parmeters
   \param data Data associated with window
   \param head Header sassoicated with window
*/
  void writeWindow(const std::vector<int> &nw, const std::vector<int> &fw,
                   const std::vector<int> &jw, const void *data,
                   void *head = 0);
  /*!
  Read a stream
   \param sz Number of bytes
   \param data Data associated with window
*/
  virtual void readStream(const long long sz, void *data) {
    throw SEPException(std::string("readStream not defined"));
  }
  /*!
Write a stream
 \param sz Number of bytes
 \param data Data associated with window
*/
  virtual void writeStream(const long long sz, const void *dat) {
    throw SEPException(std::string("writeStream not defined"));
  }
  /*!
Read a trace stream
 \param sz Number of bytes
 \param data Data associated with window
 \param head Header sossicated with stream
*/
  virtual void readTraceStream(const long long sz, void *dat, void *head = 0);
  /*!
Write a trace stream
\param sz Number of bytes
\param data Data associated with window
\param head Header sossicated with stream
*/
  virtual void writeTraceStream(const long long sz, const void *dat,
                                const void *head = 0);
  /*!
    Write reel header
    \param reelH Reel header
    */
  virtual void writeReelHead(const void *reelH);
  /*!
    Get the current position in a file
    */
  virtual long long getCurrentPos() const {
    throw SEPException(std::string("getCurrentPos not defined"));
  }
  /*!
    Get the size of the dataset
    */
  virtual long long getSize() {
    throw SEPException(std::string("getSize is undefined"));
  }
  /*!
     Swap bytes of floats
     \param n Number of floats
     \param buf Buffer to swap the bytes
     */
  void swap_float_bytes(const int n, float *buf);
  /*
    Seek to a position in a file
    \param off Offset
    \param whence from 0-begining, 1-current, 2-end of file
    */
  virtual inline void seekTo(const long long off, const int whence) {
    throw SEPException(std::string("seekTo is undefined"));
  }
  /*!
    Seek to a gicen position in a file
    \param pos Position to go to
    */
  virtual inline void seekToPos(const long long pos) {
    throw SEPException(std::string("seekToPos is undefined"));
  }
  /*!
     Set file parameters
     \param nm Name of file
     \param usage Usage for file
     \param reelH Reel header
     \param Trace header
     \param Element size of data
     \param swapData Whether or not to swap the bytes of the dataset
     \param hyper Hypercube assoicated with the dataset
     */
  void setFileParams(const std::string nm, const usage_code usage,
                     const int reelH, const int traceH, const int esize,
                     const bool swapData, std::shared_ptr<hypercube> hyper);

 private:
  /*<
    Read blocks of a dataset
      \param naxes Number of axes
      \param nwo,fwo,jwo
      */
  void readBlocks(const int naxes, const std::vector<int> &nwo,
                  const std::vector<int> &fwo, const std::vector<int> &jwo,
                  const std::vector<int> &nwi, const std::vector<int> &fwi,
                  const std::vector<int> &jwi, const long long buf, void *data,
                  void *head);
  void writeBlocks(const int naxes, const std::vector<int> &nwo,
                   const std::vector<int> &fwo, const std::vector<int> &jwo,
                   const std::vector<int> &nwi, const std::vector<int> &fwi,
                   const std::vector<int> &jwi, const long long buf,
                   const void *data, const void *head);

 protected:
  usage_code _usage;
  FILE *_myf;
  bool _swapData = false;
  long long _reelH, _traceH;
  long long _esize;
  std::string _nm;
  std::shared_ptr<hypercube> _hyper;
};
class myFileIO : public basicIO {
 public:
  myFileIO(const std::string &nm, const usage_code usage, const int reelH,
           const int traceH, const int esize, const bool swapData,
           std::shared_ptr<hypercube> hyper);
  virtual inline void seekToPos(const long long pos) {
    long long ft = getCurrentPos();
    long long bg = 1024 * 1024 * 1024;
    long long diff = pos - ft;
    while (diff != 0) {
      long long dst;
      if (diff > 0)
        dst = std::min(bg, diff);
      else
        dst = -std::min(-diff, bg);
      fseek(_myf, dst, SEEK_CUR);
      diff -= dst;
    }
  }
  virtual inline void seekTo(const long long pos, const int whence) {
    fseek(_myf, pos, whence);
  }
  virtual long long getCurrentPos() const { return ftell(_myf); }
  virtual long long getSize() {
    long long cp = getCurrentPos();
    fseek(_myf, 0L, SEEK_END);
    long long end = getCurrentPos();
    fseek(_myf, cp, 0);
    return end;
  }
  virtual void readStream(const long long sz, void *data);
  virtual void writeStream(const long long sz, const void *data);
  virtual void close() {
    if (myf != 0) fclose(myf);
  }
  ~myFileIO() { close(); }
  FILE *myf = 0;
};

void partsToBlock(const std::shared_ptr<hypercube> hyper, const int traceH,
                  const int esize, const std::vector<int> &nw,
                  const std::vector<int> &fw, const std::vector<int> &jw,
                  void *in, const void *out, const void *head);
void blockToParts(const std::shared_ptr<hypercube> hyper, const int traceH,
                  const int esize, const std::vector<int> &nw,
                  const std::vector<int> &fw, const std::vector<int> &jw,
                  const void *in, void *out, void *head);

}  // namespace SEP

#endif
