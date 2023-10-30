#ifndef DEBUG_IO_H
#define DEBUG_IO_H 1
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>  // std::stringstream
#include <string>
#include <vector>
namespace SEP {
/*!
A  virutal class for debugging. Useful for both straight C++ and python.
*/
class debug {
 public:
  debug() { ; }
  //! Add a line to buffer
  /*!
    \param line What to add to debug buffer
  */
  void addLine(std::stringstream &line);
  //! Flush the buffer
  virtual void flush() = 0;
  //! Return all the current lines in the buffer
  std::vector<std::string> returnLines() { return _lines; }
  //! Set how often to automatically flus the buffer
  void setFlush(int nflush) { _nflush = nflush; }
  //! How to close the buffer
  virtual void close() { ; }
  std::vector<std::string> _lines;  ///< Buffer of currerntly stored lines

 private:
  int _nflush;  /// How often to flush
};

/*!
Debuging class to the screen
*/
class errDebug : public debug {
 public:
  //! Create a debugging object
  errDebug(int nflush = 1000000) { setFlush(nflush); }

  virtual void flush() override;
};
/*!
Debuging class to a file
*/
class fileDebug : public debug {
 public:
  /*!
A  virutal class for debugging. Useful for both straight C++ and python.
*/
  /*!
  \param name Name of file to write info to
  \param nflush How many lines to store before automatically flushing to file
*/
  fileDebug(std::string &name, int nflush = 1000000);
  virtual void close() override;
  virtual void flush() override;

 private:
  std::ofstream _file;
};

/*!
A  singleton class to store debuggingp information
*/

class debugging {
 public:
  //! Return the instance of the debugging object

  static debugging &instance() {
    static debugging instance;
    return instance;
  }
  /*!
Add debug case
*/
  /*!
  \param name Name of the debugging info
  \param deb Debugging class
*/
  void addDebug(std::string name, std::shared_ptr<debug> deb) {
    _debugs[name] = deb;
  }

  /*!
Get debugging info
*/
  /*!
  \param name Name of the debugging info
*/
  std::shared_ptr<debug> getDebug(const std::string &name);

 private:
  debugging() { ; }
  ~debugging() { ; }
  std::map<std::string, std::shared_ptr<debug>> _debugs;  /// Debugging info
};

}  // namespace SEP
   /*!
Return debug buffer
*/

std::vector<std::string> getDebug(std::string buffer);
#endif
