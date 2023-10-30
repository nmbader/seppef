#ifndef IRREG_FILE_H
#define IRREG_FILE_H
#include "header.h"
#include "hypercube.h"
/*
  Regular data is n-d hypercube (for regular seismic a single trace)
    Each regular n-d cube has a header associated with it
    Headers exist in a grid (could be as simple at 1-D list (collection of
  traces)) Headers potentially in blocks. Each block knows how to access the
  traces assoicated with Grid knows how to access header blocks Traces can be in
  order or out of order. In order more possibilites to do combining operations

    When you do a read of headers you are returned an object *headerBlock* that
  knows how to read regular cubes associated with it.

    You can then specify to read some subset of the headers, which will return
  you data associated with it

    You create a headerBlock by passing the headers and traces
   */

namespace SEP {

class IOHeaderBlock {
 public:
  IOHeaderBlock()
};

class irregFile {
 public:
 private:
  std::shared_ptr<hypercube> _irregHyper;
  std::shared_ptr<hypercube> _regHyper;
  std::vector<key> _keys;
};

}  // namespace SEP
#endif