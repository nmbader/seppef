#ifndef vector_h
#define vector_h 1
#include <cstdint>
#include <memory>
#include <sstream>
#include "SEPException.h"
namespace SEP {
/*!
 Generic vector class. Useful for inversion and other applications.
*/
class Vector {
 public:
  //! Create a vector (virtual class)
  Vector() : _spaceOnly(false), _vecType(std::string("base")) { ; }
  //! Add vector to another vector  self+=vec
  /*!
    \param vec Vector to add to the current vector
  */
  virtual void add(const std::shared_ptr<Vector> vec) {
    throw SEPException(std::string("Must override add"));
  }
  //! Scale vector self*=scale
  /*!
    \param scale What to scale vector by
  */

  virtual void scale(const double scale) {
    throw SEPException(std::string("Must override scale"));
  }

  //! Add to vector scaling each self=self*sc1+vec*sc2
  /*!
    \param sc1 What to scale current vector by
    \param vec Other vector to scale add to current vector
    \param sc2 What to scale the second vector by
  */

  virtual void scaleAdd(const double sc1, const std::shared_ptr<Vector> vec,
                        const double sc2) {
    throw SEPException(std::string("Must override scaleAdd"));
  }
  //! Fill vector with random number
  virtual void random() {
    throw SEPException(std::string("Must override random"));
  }
  /*! Apply a softclip to vector
       vec=sc*vec /sqrt(1+sc*sc*vec*vec)

       \param sc Value to use in softclip
  */
  virtual void softClip(float sc) {
    throw SEPException(std::string("Must override softclip"));
  }
  //! Return the absolute maximum value of vector
  virtual double absMax() const {
    throw SEPException(std::string("Must override absMax"));
  }
  /*!  Return the dot product of current vec and another vec
     return SUM(this[]*vec2[])
     \param vec2 Vector to calculate the dot prodcut with
     */
  virtual double dot(const std::shared_ptr<Vector> vec2) const {
    throw SEPException(std::string("Must override dot"));
  }

  /*!  Check to see if current vector belongs to the same space as vec2

   \param vec2 Vector to check the space with
   \param checkAlloc Also check that the vector is allocated
   */
  virtual bool checkSame(const std::shared_ptr<Vector> vec2,
                         const bool checkAlloc = false) const {
    throw SEPException(std::string("Must override checkSame"));
  }
  /*! Remove storage for vector, retain only vector space*/
  virtual void cleanMemory() {
    throw SEPException(std::string("Must override cleanMemory"));
  }
  /*! Vector only contains space representation*/
  virtual void setSpace() { _spaceOnly = true; }

  /*!  Return  information about vector (debugging)

   \param lev  Level of debugging information to provide
   \param str  Stream to add debugging info to
   */
  virtual void infoStream(const int lev, std::stringstream &str) {
    throw SEPException(std::string("Must override infoStream"));
  }
  /*!  Return  information about vector (debugging)

  \param nm   Mark this debugging information nm
  \param lev  Level of debugging information to provide
  */
  std::string info(const std::string &nm, const int lev);

  //! Calculate checksum for data
  virtual void calcCheckSum() {
    throw SEPException(std::string("Must override calcCheckSum"));
  }
  /*!  Store checksum va;iue

  \param val Value of checksum to store
  */
  void setCheckSum(const uint64_t val) { _checkSum = val; }
  ///! Return checksum value
  uint64_t getCheckSum() { return _checkSum; }
  ///! Return whether this vector only represents space
  inline bool spaceOnly() const { return _spaceOnly; }
  ///! Set that this vector contains storage for vector elements
  void setNotSpace() { _spaceOnly = false; }
  std::string _vecType;  ///< Vector type

 protected:
  bool _spaceOnly;     ///< Whether or not vector is only  contains vector space
  uint64_t _checkSum;  ///< Checksum for vector
};

}  // namespace SEP
#endif
