#ifndef int1d_reg_h
#define int1d_reg_h 1
#include <intHyper.h>
#include <cstdint>
#include <iostream>
#include "boost/multi_array.hpp"
#include "int2DReg.h"
#include "int3DReg.h"
#include "int4DReg.h"
#include "int5DReg.h"
#include "int6DReg.h"

typedef boost::multi_array<int, 1> int1D;
namespace SEP {
/*!
A regular sampled 1-D function with integer storage
*/
class int1DReg : public intHyper {
 public:
  /*!
 Create a 1-D integer vector from a hypercube
      \param hyper Hypercube describing RSF

 */
  int1DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 1-D integer vector from just lengths
    \param n  Dimensions of the hypercube

*/
  int1DReg(const int n) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
  Create a 1-D integer vector from axes
       \param a1 Axes of the hypercube

  */
  int1DReg(const SEP::axis &a1) {
    std::vector<SEP::axis> as(1, a1);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    initNoData(hyp);
  }
  /*!
 Create a 1-D integer vector from just lengths
      \param hyper Hypercube for RSF
     \param vals Values to fill vector with

 */
  int1DReg(std::shared_ptr<SEP::hypercube> hyper, const int1D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 1-D integer vector from just lengths
  \param n  Dimensions of the hypercube
 \param vals Values to fill vector with

*/
  int1DReg(const int n, int1D &vals) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
Create a 1-D integer vector from axes
     \param a1 Axes of the hypercube
    \param vals Values to fill vector with

*/
  int1DReg(const SEP::axis &a1, const int1D &vals) {
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
  int1DReg(const std::shared_ptr<int6DReg> old, const int iax1, const bool rev1,
           const std::vector<int> &ipos, const std::vector<int> &beg,
           const std::vector<int> &end);
  /*!
    Create 1-D byte vector from a subset of byte RSF

    \param old RSF to subsample
    \param iax1 Fast axis of the output
    \param rev1 Whether or not to reverse the axis
    \param ipos Location if the subset to grab
    \param beg,end Begining and end of the subset to grab along
    */
  int1DReg(const std::shared_ptr<int5DReg> old, const int iax1, const bool rev1,
           const std::vector<int> &ipos, const std::vector<int> &beg,
           const std::vector<int> &end);
  /*!
   Create 1-D byte vector from a subset of byte RSF

   \param old RSF to subsample
   \param iax1 Fast axis of the output
   \param rev1 Whether or not to reverse the axis
   \param ipos Location if the subset to grab
   \param beg,end Begining and end of the subset to grab along
   */
  int1DReg(const std::shared_ptr<int4DReg> old, const int iax1, const bool rev1,
           const std::vector<int> &ipos, const std::vector<int> &beg,
           const std::vector<int> &end);
  /*!
   Create 1-D byte vector from a subset of byte RSF

   \param old RSF to subsample
   \param iax1 Fast axis of the output
   \param rev1 Whether or not to reverse the axis
   \param ipos Location if the subset to grab
   \param beg,end Begining and end of the subset to grab along
   */
  int1DReg(const std::shared_ptr<int3DReg> old, const int iax1, const bool rev1,
           const std::vector<int> &ipos, const std::vector<int> &beg,
           const std::vector<int> &end);
  /*!
   Create 1-D byte vector from a subset of byte RSF

   \param old RSF to subsample
   \param iax1 Fast axis of the output
   \param rev1 Whether or not to reverse the axis
   \param ipos Location if the subset to grab
   \param beg,end Begining and end of the subset to grab along
   */
  int1DReg(const std::shared_ptr<int2DReg> old, const int iax1, const bool rev1,
           const std::vector<int> &ipos, const std::vector<int> &beg,
           const std::vector<int> &end);
  /*!
   Create 1-D byte vector from a subset of byte RSF

   \param old RSF to subsample
   \param iax1 Fast axis of the output
   \param rev1 Whether or not to reverse the axis
   \param ipos Location if the subset to grab
   \param beg,end Begining and end of the subset to grab along
   */
  int1DReg(const std::shared_ptr<int1DReg> old, const int iax1, const bool rev1,
           const std::vector<int> &ipos, const std::vector<int> &beg,
           const std::vector<int> &end);
  /*!
   Allocate data for vector
   */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new int1D(boost::extents[ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
     Return a subset of the vector
     \param nw,fw,jw Windowing parameters
     */
  std::shared_ptr<int1DReg> window(const std::vector<int> &nw,
                                   const std::vector<int> &fw,
                                   const std::vector<int> &jw) const;
  /*!
    Return a subset of the vector
    \param nw,fw,jw Windowing parameters
    */
  std::shared_ptr<int1DReg> window(const int nw, const int fw, const int jw) {
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
  std::shared_ptr<int1DReg> clone() const;
  /*!
  Make a copy of the vector space
  */
  std::shared_ptr<int1DReg> cloneSpace() const;
  /*!
Deallocate storage for vector, turn into vector space
*/
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<int1D> _mat;  ///< Boost storage for vector

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

  void initData(std::shared_ptr<SEP::hypercube> hyper, const int1D &vals);
};

}  // namespace SEP
#endif
