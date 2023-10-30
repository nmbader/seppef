#ifndef reg_space_h
#define reg_space_h 1
#include "SEPException.h"
#include "hypercube.h"

namespace SEP {
/*!
  A class for regular sampled function (virtual)
*/
class regSpace {
 public:
  /*!
     Initalizer for virtual regular sampled function class
  */
  regSpace() { ; }
  /*!
     Return the hypercube associated with the regular sampled function
*/
  std::shared_ptr<SEP::hypercube> getHyper() const { return _hyper; }
  /*!
     Set the hypercube that describes the regular sampled function
    */
  void setHyper(std::shared_ptr<SEP::hypercube> h) { _hyper = h->clone(); }
  /*!
      Set the memory associted with regular sampled function and the element
     size of the data

      \param ptr  Pointer associated with the data
      \param esize Element size sizeof(float),sizeof(int),etc
  */
  void setMemPtr(void *ptr, size_t esize) {
    _storage = ptr;
    _esize = esize;
  }
  /*!
     Checck to see if window parameters make sense

     \param n  Axis length
     \param nw Number of samples to grab
     \param fw First element to grab
     \param jw Samples to skip
   */

  void checkWindow(const int n, const int nw, const int fw, const int jw) const;
  /*!
    Given hypercube translate an axis to a key

    \param iax Axis to convert to a key
    */
  std::vector<float> axisToKey(const int iaxis) const;

  ///*  Return pointer to storage
  void *getVoidPtr() { return _storage; }
  ///* Get the size of each element in  vector
  size_t getEsize() { return _esize; }
  ///* Remove the regular space object
  virtual ~regSpace() = default;

 protected:
  /*!
      Return how to loop through a dataset to return a slice with arbitray
     orientations

      \param  nd of the dataset
      \param iax1 Fast axis of output traverse
      \param rev1 Whether or not reverse fast axis
      \param f1  First element to grab along fast axis
      \param j1  Skip every j1 elements along fast axis
      \param iax2 Slow axis of output trarvers
      \param rev2 Whether or not to reverse slow axis
      \param f2   First element of slow axis
      \param j2 Skip every j2 element along slow axis
  */
  inline void calcTraverse(const std::vector<int> &nd, const int iax1,
                           const bool rev1, int &f1, int &j1, const int iax2,
                           const bool rev2, int &f2, int &j2) {
    if (iax2 > iax1) {
      if (!rev1) {
        f1 = 0;
        j1 = 1;
      } else {
        f1 = nd[iax1] - 1;
        j1 = -1;
      }
      if (!rev2) {
        f2 = 0;
        j2 = nd[iax1];
      } else {
        f2 = (nd[iax2] - 1) * nd[iax1];
        j2 = -nd[iax1];
      }
    } else {
      if (!rev2) {
        f2 = 0;
        j2 = 1;
      } else {
        f2 = nd[iax2] - 1;
        j2 = -1;
      }
      if (!rev1) {
        f1 = 0;
        j1 = nd[iax2];
      } else {
        f1 = (nd[iax1] - 1) * nd[iax2];
        j1 = -nd[iax2];
      }
    }
  }

 private:
  std::shared_ptr<SEP::hypercube>
      _hyper;      ///< Hypercube associated with the regular sampled function
  size_t _esize;   ///< Element size for each vector element
  void *_storage;  ///< Pointer to storage for regular sampled function
};
}  // namespace SEP
#endif
