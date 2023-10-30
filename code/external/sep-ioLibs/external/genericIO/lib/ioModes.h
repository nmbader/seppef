#ifndef IOMODES_H
#define IOMODES_H 1
#include <string.h>
#include <mutex>
#include "genericIO.h"
#include "jsonGenericIO.h"
namespace SEP {
/*!
  Class controlling different IO methods

  */
class ioModes {
 public:
  /*!
    Initialize default IO mode
    */
  ioModes() { ; }
  /*!
    Initialize IOmodes from a vector or arguments
    \param args Command line arguments
    */
  ioModes(std::vector<std::string> args) {
    int nc = args.size();
    char **argv;
    if (nc == 0) {
      nc = 1;
      argv = new char *[1];
      argv[0] = new char[1];
      strcpy(argv[0], "");

    } else {
      argv = new char *[args.size()];
      for (auto i = 0; i < args.size(); i++) {
        argv[i] = new char[args[i].size() + 1];
        strcpy(argv[i], args[i].c_str());
      }
    }
    setup(nc, argv);
    /*
    for (auto i = 0; i < args.size(); i++) {
      delete[] argv[i];
    }
    delete[] argv;
    */
  }
  /*!
     Initalize IOModes from standard C arguments
     \param argc Number of arguments
     \param argv List of arguments
     */
  ioModes(const int argc, char **argv) { setup(argc, argv); }
  /*!
     Setup IO modes arguments
     */
  void setup(const int argc, char **argv);
  /*!
    Get the default IO type
  */
  std::shared_ptr<genericIO> getDefaultIO();
  /*!
    Get a specific IO
    \param def IO type to grab
    */
  std::shared_ptr<genericIO> getIO(const std::string &def);

  /*!
    Get all of the IOs available

    @return list of IOs
  */
  std::vector<std::string> getIOs() const;
  /*!
     Return a file object
     \param tag Tag for file
     \param ioname Name of IO
     \param name Name of file
     \param usage Usage for  file
     */

  std::shared_ptr<genericRegFile> getRegFileTag(const std::string &tag,
                                                const std::string &ioname,
                                                const std::string &name,
                                                const usage_code usage);
  /*!
         Return IO user wants to use for input
  */
  std::shared_ptr<genericIO> getInputIO();

  /*!
         Return IO user wants to use for output
  */
  std::shared_ptr<genericIO> getOutputIO();

  /*!
         Return paramObj  user wants to use for parameter handle
  */
  std::shared_ptr<paramObj> getParamObj();

  /*!
     Set param object used for all IO objects
   */
  void changeParamObj(std::shared_ptr<paramObj> par);

  /*!
    Return file object from default IO
    \param name Name of file
    \param Usage for file
    */

  std::shared_ptr<genericRegFile> getGenericRegFile(const std::string &name,
                                                    const usage_code usage);
  /*
    Return the default IO type
    */
  std::string getDefaultType() { return _defaultType; }

 private:
  std::shared_ptr<genericIO> _defaultIO;  ///< Pointer to the default IO
  std::map<std::string, std::shared_ptr<genericIO>>
      _ios;  ///< IO type dictionary

  std::shared_ptr<paramObj> _par;  ///< Parameter object
  std::string _defaultType;        ///< Default type
};
/*!
  Class for accessing IO from fortran. Singleton class
  */
class ioModesFortran {
 private:
  /*!
    Initialize default fortran object
    */
  ioModesFortran() { ; }

  ioModesFortran(const ioModesFortran &rs);
  ioModesFortran &operator=(const ioModesFortran &rs);
  static std::shared_ptr<ioModesFortran> instance;

  std::shared_ptr<ioModes> _io;

 public:
  void setup(const int argc, char **argv);
  ~ioModesFortran() { ; }
  std::shared_ptr<ioModes> getIOModes() { return _io; }

  static std::shared_ptr<ioModesFortran> &getInstance() {
    if (!instance) {
      instance.reset(new ioModesFortran());
      instance->_io.reset(new ioModes());
    }

    return instance;
  }
};
}  // namespace SEP

#endif
