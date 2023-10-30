#ifndef int2d_reg_h
#define int2d_reg_h 1
#include <intHyper.h>
#include "int3DReg.h"
#include "int4DReg.h"
#include "int5DReg.h"
#include "int6DReg.h"

#include "boost/multi_array.hpp"

typedef boost::multi_array<int, 2> int2D;
namespace SEP {
/*!
A regular sampled 2-D function with integer storage
*/
class int2DReg : public intHyper {
 public:
  /*!
   Create a 2-D integer vector from a hypercube
        \param hyper Hypercube describing RSF

   */
  int2DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 2-D integer vector from just lengths
     \param n1,n2 Dimensions of the hypercube

*/
  int2DReg(const int n1, const int n2) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 2-D integer vector from axes
      \param a1,a2 Axes if the hypercube

 */
  int2DReg(const SEP::axis &a1, const SEP::axis &a2) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 2-D integer vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Data to fill vector with

 */
  int2DReg(std::shared_ptr<SEP::hypercube> hyper, const int2D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 2-D integer vector from just lengths
  \param n1,n2 Dimensions of the hypercube
  \param vals Data to fill vector with

*/
  int2DReg(int n1, int n2, int2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
 Create a 2-D integer vector from axes
      \param a1,a2 Axes if the hypercube
      \param vals Data to fill vector with
 */
  int2DReg(SEP::axis &a1, SEP::axis &a2, int2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
 Create 2-D byte vector from a subset of byte RSF

 \param old RSF to subsample
 \param iax1 Fast axis of the output
 \param rev1 Whether or not to reverse the  first axis
 \param iax2 Slower axis of the output
 \param rev2 Whether or not reverse the second axis
 \param ipos Location if the subset to grab
 \param beg,end Begining and end of the subset to grab along
 */
  int2DReg(const std::shared_ptr<int2DReg> old, const int iax1, const bool rev1,
           const int iax2, const bool rev2, const std::vector<int> &ipos,
           const std::vector<int> &beg, std::vector<int> &end);
  /*!
  Create 2-D byte vector from a subset of byte RSF

  \param old RSF to subsample
  \param iax1 Fast axis of the output
  \param rev1 Whether or not to reverse the  first axis
  \param iax2 Slower axis of the output
  \param rev2 Whether or not reverse the second axis
  \param ipos Location if the subset to grab
  \param beg,end Begining and end of the subset to grab along
  */
  int2DReg(const std::shared_ptr<int3DReg> old, const int iax1, const bool rev1,
           const int iax2, const bool rev2, const std::vector<int> &ipos,
           const std::vector<int> &beg, std::vector<int> &end);
  /*!
  Create 2-D byte vector from a subset of byte RSF

  \param old RSF to subsample
  \param iax1 Fast axis of the output
  \param rev1 Whether or not to reverse the  first axis
  \param iax2 Slower axis of the output
  \param rev2 Whether or not reverse the second axis
  \param ipos Location if the subset to grab
  \param beg,end Begining and end of the subset to grab along
  */
  int2DReg(const std::shared_ptr<int4DReg> old, const int iax1, const bool rev1,
           const int iax2, const bool rev2, const std::vector<int> &ipos,
           const std::vector<int> &beg, std::vector<int> &end);
  /*!
 Create 2-D byte vector from a subset of byte RSF

 \param old RSF to subsample
 \param iax1 Fast axis of the output
 \param rev1 Whether or not to reverse the  first axis
 \param iax2 Slower axis of the output
 \param rev2 Whether or not reverse the second axis
 \param ipos Location if the subset to grab
 \param beg,end Begining and end of the subset to grab along
 */
  int2DReg(const std::shared_ptr<int5DReg> old, const int iax1, const bool rev1,
           const int iax2, const bool rev2, const std::vector<int> &ipos,
           const std::vector<int> &beg, std::vector<int> &end);
  /*!
 Create 2-D byte vector from a subset of byte RSF

 \param old RSF to subsample
 \param iax1 Fast axis of the output
 \param rev1 Whether or not to reverse the  first axis
 \param iax2 Slower axis of the output
 \param rev2 Whether or not reverse the second axis
 \param ipos Location if the subset to grab
 \param beg,end Begining and end of the subset to grab along
 */
  int2DReg(const std::shared_ptr<int6DReg> old, const int iax1, const bool rev1,
           const int iax2, const bool rev2, const std::vector<int> &ipos,
           const std::vector<int> &beg, std::vector<int> &end);
  /*!
  Allocate data for vector
  */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new int2D(boost::extents[ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
   Make a copy of the vector
   */
  std::shared_ptr<int2DReg> clone() const;
  /*!
 Make a copy of the vector space
 */
  std::shared_ptr<int2DReg> cloneSpace() const;
  /*!
  Create a window of the current vector
  \param nw,fw,jw Window parameters
  */
  std::shared_ptr<int2DReg> window(const std::vector<int> &nw,
                                   const std::vector<int> &fw,
                                   const std::vector<int> &jw) const;
  /*!
 Deallocate storage for vector, turn into vector space
  */
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<int2D> _mat;  ///< Storage for vector

 private:
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
  void initData(std::shared_ptr<SEP::hypercube> hyper, const int2D &vals);
};
}  // namespace SEP
#endif
