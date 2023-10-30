#ifndef buffers_irreg_file_h
#define buffers_irreg_file_h 1
#include <stdbool.h>
#include <string>
#include "basicIO.h"
#include "buffers.h"
#include "json.h"
#include "jsonGenericIrregFile.h"
namespace SEP {

/*!
  Object for handling irregular files broken into blocks
*/
class buffersIrregFile : public jsonGenericIrregFile {
 public:
  /*!
     Initialize an irregular file
 */
  buffersIrregFile() { ; }
  /*!
     Write the description for the file
 */
  virtual void writeDescription();
  /*!
     Set how memory should be handled for buffer object
       \param mem Memory usage object
 */
  void setMemoryUsage(std::shared_ptr<SEP::IO::memoryUsage> mem) {
    if (!_hyper) throw SEPException(std::string("Hypercube has not been set"));
    _mem = mem;
  }
  /*!
   Set the compression scheme for the dataset
     \param com Compresison object
*/
  void setCompression(std::shared_ptr<SEP::IO::compress> com) { _comp = com; }
  /*!
   Set how the dataset should be broken into blocks
     \param block Blocking object
*/
  void setBlocking(std::shared_ptr<SEP::IO::blocking> block) { _block = block; }
  /*!
     Create the buffers
 */
  virtual void createBuffers() = 0;

 protected:
  std::shared_ptr<SEP::IO::buffers> _bufs = nullptr;  ///< Buffering scheme
  std::shared_ptr<SEP::IO::memoryUsage> _mem =
      nullptr;                                         ///< Memory usage scheme
  std::shared_ptr<SEP::IO::compress> _comp = nullptr;  ///< Compression scheme
  std::shared_ptr<SEP::IO::blocking> _block = nullptr;

};  // namespace SEP

}  // namespace SEP
#endif
