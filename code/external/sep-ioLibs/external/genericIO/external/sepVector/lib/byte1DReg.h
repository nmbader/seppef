#ifndef byte1d_reg_h
#define byte1d_reg_h 1
#include <byteHyper.h>
#include <cstdint>
#include <iostream>
#include "boost/multi_array.hpp"
#include "byte2DReg.h"
#include "byte3DReg.h"
#include "byte4DReg.h"
#include "byte5DReg.h"
#include "byte6DReg.h"

typedef boost::multi_array<unsigned char, 1> byte1D;
namespace SEP {
/*!
A regular sampled 1-D function with unsigned char storage
*/
class byte1DReg : public byteHyper {
 public:
  /*!
   Create a 1-D unsigned char vector from a hypercube
        \param hyper describing RSF

   */
  byte1DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
 Create a 1-D unsigned char vector from just lengths
      \param n  Dimensions of the hypercube

 */

  byte1DReg(const int n) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
   Create a 1-D unsigned char vector from axes
        \param a1 Axes of the hypercube

   */
  byte1DReg(const SEP::axis &a1) {
    std::vector<SEP::axis> as(1, a1);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    initNoData(hyp);
  }
  /*!
   Create a 1-D unsigned char vector from just lengths
        \param hyper  Hypercube describing RSF
       \param vals Values to fill vector with

   */
  byte1DReg(std::shared_ptr<SEP::hypercube> hyper, const byte1D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 1-D unsigned char vector from just lengths
    \param n  Dimensions of the hypercube
   \param vals Values to fill vector with

*/
  byte1DReg(const int n, byte1D &vals) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
  Create a 1-D unsigned char vector from axes
       \param ax Axes of the hypercube
      \param vals Values to fill vector with

  */
  byte1DReg(const SEP::axis &ax, const byte1D &vals) {
    std::vector<SEP::axis> as(1, ax);
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
  byte1DReg(const std::shared_ptr<byte6DReg> old, const int iax1,
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
  byte1DReg(const std::shared_ptr<byte5DReg> old, const int iax1,
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
  byte1DReg(const std::shared_ptr<byte4DReg> old, const int iax1,
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
  byte1DReg(const std::shared_ptr<byte3DReg> old, const int iax1,
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
  byte1DReg(const std::shared_ptr<byte2DReg> old, const int iax1,
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
  byte1DReg(const std::shared_ptr<byte1DReg> old, const int iax1,
            const bool rev1, const std::vector<int> &ipos,
            const std::vector<int> &beg, const std::vector<int> &end);
  /*!
    Allocate data for vector
    */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new byte1D(boost::extents[ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
  Make a copy of the vector
  */
  std::shared_ptr<byte1DReg> clone() const;
  /*!
  Make a copy of the vector
  */
  std::shared_ptr<byte1DReg> cloneSpace() const;
  /*!
 Deallocate storage for vector, turn into vector space
  */
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<byte1D> _mat;  ///< Storage for vector
  /*!
     Return a subset of the vector
     \param nw,fw,jw Windowing parameters
     */
  std::shared_ptr<byte1DReg> window(const std::vector<int> &nw,
                                    const std::vector<int> &fw,
                                    const std::vector<int> &jw) const;
  /*!
    Return a subset of the vector
    \param nw,fw,jw Windowing parameters
    */
  std::shared_ptr<byte1DReg> window(const int nw, const int fw, const int jw) {
    std::vector<int> nws;
    nws.push_back(nw);
    std::vector<int> fws;
    fws.push_back(fw);
    std::vector<int> jws;
    jws.push_back(jw);
    return (window(nws, fws, jws));
  }

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
  void initData(std::shared_ptr<SEP::hypercube> hyper, const byte1D &vals);
};

}  // namespace SEP
#endif
