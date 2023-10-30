#ifndef byte_hyper_h
#define byte_hyper_h 1
#include <hypercube.h>
#include <cstdint>
#include <sstream>
#include "int1DReg.h"
#include "regSpace.h"
namespace SEP {
/*!
  A regular sampled function that stores unsigned char values. Storage is
  actually done in inherited classes (1-D,2-D,3-D)
*/
class byteHyper : public regSpace {
 public:
  /*!
Initializer for byteHyper class. Only used by inherited class
*/
  byteHyper() { ; }
  /*!
Initializer for byteHyper class. Inititalize form a hypercube
  \param hyper Hypercube describing the RSF
*/

  byteHyper(std::shared_ptr<SEP::hypercube> hyper) { setHyper(hyper->clone()); }
  /*!
     Return whether or not this just containing space information
  */
  //! Fill vector with random number
  virtual void random();

  /*
    \param mn  Minimum value for histogram
    \param mx  Maximum value for histogram
    \param histo Number of elements in the histogram
   */
  virtual void calcHisto(std::shared_ptr<int1DReg> &histo, float mn, float mx);
  /*!  Set that this is only vector space with no storage
   */
  virtual void setSpace() { _spaceOnly = true; }
  /*!  Set that this has storage along with the vector sapce
   */
  virtual void setNotSpace() { _spaceOnly = false; }
  /*!
   Return whether or not this just containing space information
*/
  inline bool getSpaceOnly() const { return _spaceOnly; }
  /*!  Return the pct value of the dataset

    \param pct Percentage value to return
    \param j   Approximate the answer by ony taking every jth value

   */
  unsigned char cent(const float pct, const int j) const {
    long long iv = std::max(
        (long long)0, std::min((long long)(getHyper()->getN123() * pct / 100.),
                               getHyper()->getN123() - 1));
    return cent(iv, j);
  }
  /*!  Return the iv value of the dataset if sorted low to high

    \param iv  Value of sorted data to return
    \param j   Approximate the answer by ony taking every jth value

   */
  unsigned char cent(const long long iv, const int j) const;
  /*!
    Clip a dataset
        a[]=std::min(eclip,max(bclip,a[]))

        or

        a[]=std::max(eclip,min(bclip,a[]))

    \param  bclip Minimum value
    \param  eclip Maximum value

*/
  void clip(const unsigned char bclip, const unsigned char eclip);
  /*!
   Set a pointer to the storage for the dataset.
   Storage is handled by child classes.

   \param ptr  Pointer to allocated memory
*/
  void setData(unsigned char *ptr) {
    _vals = ptr;
    setNotSpace();
    setMemPtr((void *)ptr, sizeof(unsigned char));
  }
  //! Calculate checksum for data

  void calcCheckSum();
  /*!  Store checksum va;iue

\param val Value of checksum to store
*/
  void setCheckSum(const uint64_t val) { _checkSum = val; }
  /*!
   Whether or not the current vector exists in a different space

   \param vec2 Vector space to compare to
*/
  bool isDifferent(std::shared_ptr<byteHyper> vec2) {
    calcCheckSum();
    if (vec2->getCheckSum() != getCheckSum()) return true;
    return false;
  }
  /*!
     Return the norm of the dataset.
       1 - sum(|a[]|)
       2 - 1/2(a[]*a[])

     \param nrm Norm to calculate
  */
  long long norm(const int nrm) const;
  /*!
Set the valule of vector to 0

*/
  void zero() { set(0); }
  /*!
   Set the valule of vector to a given value

   \param val Value to set vector to
 */
  void set(const unsigned char val);
  /*!
      Return the pointer to the memory for the vector
  */
  unsigned char *getVals() { return _vals; }
  /*!
  Return the pointer to the memory for the vector with const tag
*/
  const unsigned char *getCVals() const { return _vals; }
  //! Return the absolute maximum value of vector

  virtual int absMax() const;
  //! Return the  minimum value of vector

  int min() const;
  //! Return the  maximum value of vector

  int max() const;
  /*!  Return  information about vector (debugging)

\param lev  Level of debugging information to provide
\param str  Stream to add debugging info to
*/
  virtual void infoStream(const int lev, std::stringstream &str);
  /*!  Check to see if current vector belongs to the same space as vec2

   \param vec2 Vector to check the space with
   */

  virtual bool checkSame(const std::shared_ptr<SEP::byteHyper> vec2) const;
  ///! Return checksum value

  uint64_t getCheckSum() { return _checkSum; }

  bool _spaceOnly;       ///< Whether or not vector is only vector space
  std::string _vecType;  ///< Name for vector type

 private:
  unsigned char *_vals;  ///< Storage for vector
  uint64_t _checkSum;    ///< Checksum value
};

}  // namespace SEP
#endif
