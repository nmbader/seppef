#ifndef complex1d_reg_h
#define complex1d_reg_h 1
#include <complexHyper.h>
#include <cstdint>
#include <iostream>
#include "boost/multi_array.hpp"
#include "complex2DReg.h"
#include "complex3DReg.h"
#include "complex4DReg.h"
#include "complex5DReg.h"
#include "complex6DReg.h"

typedef boost::multi_array<std::complex<float>, 1> complex1D;
namespace SEP {
/*!
A regular sampled 1-D function with complex float storage
*/
class complex1DReg : public complexHyper {
 public:
  /*!
 Create a 1-D complex float vector from a hypercube
      \param hyper Hypercube describing RSF

 */
  complex1DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 1-D complex float vector from just lengths
    \param n  Dimensions of the hypercube

*/

  complex1DReg(const int n) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 1-D complex float vector from axes
      \param a1 Axes of the hypercube

 */
  complex1DReg(const SEP::axis &a1) {
    std::vector<SEP::axis> as(1, a1);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    initNoData(hyp);
  }
  /*!
 Create a 1-D complex float vector from just lengths
      \param hyper Hypercube describing RSF
     \param vals Values to fill vector with

 */
  complex1DReg(std::shared_ptr<SEP::hypercube> hyper, const complex1D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 1-D complex float vector from just lengths
  \param n  Dimensions of the hypercube
 \param vals Values to fill vector with

*/
  complex1DReg(const int n, complex1D &vals) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
Create a 1-D complex float vector from axes
     \param a1 Axes of the hypercube
    \param vals Values to fill vector with

*/
  complex1DReg(const SEP::axis &a1, const complex1D &vals) {
    std::vector<SEP::axis> as(1, a1);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    initData(hyp, vals);
  }
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  complex1DReg(const std::shared_ptr<complex6DReg> old, const int iax1,
               const bool rev1, const std::vector<int> &ipos,
               const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  complex1DReg(const std::shared_ptr<complex5DReg> old, const int iax1,
               const bool rev1, const std::vector<int> &ipos,
               const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  complex1DReg(const std::shared_ptr<complex4DReg> old, const int iax1,
               const bool rev1, const std::vector<int> &ipos,
               const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  complex1DReg(const std::shared_ptr<complex3DReg> old, const int iax1,
               const bool rev1, const std::vector<int> &ipos,
               const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  complex1DReg(const std::shared_ptr<complex2DReg> old, const int iax1,
               const bool rev1, const std::vector<int> &ipos,
               const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  complex1DReg(const std::shared_ptr<complex1DReg> old, const int iax1,
               const bool rev1, const std::vector<int> &ipos,
               const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Allocate data for vector
    */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new complex1D(boost::extents[ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
   Return a subset of the vector
   \param nw,fw,jw Windowing parameters
   */
  std::shared_ptr<complex1DReg> window(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw) const;
  /*!
   Return a subset of the vector
   \param nw,fw,jw Windowing parameters
   */
  std::shared_ptr<complex1DReg> window(const int nw, const int fw,
                                       const int jw) {
    std::vector<int> nws;
    nws.push_back(nw);
    std::vector<int> fws;
    fws.push_back(fw);
    std::vector<int> jws;
    jws.push_back(jw);
    return (window(nws, fws, jws));
  }
  /*!
 Make a copy of the vector
 */
  std::shared_ptr<complex1DReg> clone() const;
  /*!
 Make a copy of the vector space
 */
  std::shared_ptr<complex1DReg> cloneSpace() const;
  /*!
Deallocate storage for vector, turn into vector space
 */
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<complex1D> _mat;  ///< Boost storage object

 protected:
  /*!
  Initialize without data
  \param hyper Hypercube describing space
*/
  void initNoData(std::shared_ptr<SEP::hypercube> hyper);
  /*!
    Initialize with data
    \param hyper Hypercube describing space
    \param vals Data to copy in
 */
  void initData(std::shared_ptr<SEP::hypercube> hyper, const complex1D &vals);
};

}  // namespace SEP
#endif
