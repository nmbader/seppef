#ifndef float2d_reg_h
#define float2d_reg_h 1
#include <float3DReg.h>
#include <float4DReg.h>
#include <float5DReg.h>
#include <float6DReg.h>
#include <floatHyper.h>

#include "boost/multi_array.hpp"
typedef boost::multi_array<float, 2> float2D;
namespace SEP {
/*!
A regular sampled 2-D function with float storage
*/
class float2DReg : public floatHyper {
 public:
  /*!
    Default object
    */
  float2DReg() { ; }

  /*!
      Create a 2-D float vector from a hypercube
           \param hyper Hypercube describing RSF

      */
  float2DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 2-D float vector from just lengths
     \param n1,n2 Dimensions of the hypercube

*/
  float2DReg(const int n1, const int n2) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
  Create a 2-D float vector from axes
       \param a1,a2 Axes if the hypercube

  */
  float2DReg(const SEP::axis &a1, const SEP::axis &a2) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 2-D float vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Data to fill vector with

 */
  float2DReg(std::shared_ptr<SEP::hypercube> hyper, const float2D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 2-D float vector from just lengths
  \param n1,n2 Dimensions of the hypercube
  \param vals Data to fill vector with

*/
  float2DReg(int n1, int n2, float2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
 Create a 2-D float vector from axes
      \param a1,a2 Axes if the hypercube
      \param vals Data to fill vector with
 */
  float2DReg(SEP::axis &a1, SEP::axis &a2, float2D &vals) {
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
  float2DReg(const std::shared_ptr<float6DReg> old, const int iax1,
             const bool rev1, const int iax2, const bool rev2,
             const std::vector<int> &ipos, const std::vector<int> &beg,
             const std::vector<int> &end);
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
  float2DReg(const std::shared_ptr<float5DReg> old, const int iax1,
             const bool rev1, const int iax2, const bool rev2,
             const std::vector<int> &ipos, const std::vector<int> &beg,
             const std::vector<int> &end);
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
  float2DReg(const std::shared_ptr<float4DReg> old, const int iax1,
             const bool rev1, const int iax2, const bool rev2,
             const std::vector<int> &ipos, const std::vector<int> &beg,
             const std::vector<int> &end);
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
  float2DReg(const std::shared_ptr<float3DReg> old, const int iax1,
             const bool rev1, const int iax2, const bool rev2,
             const std::vector<int> &ipos, const std::vector<int> &beg,
             const std::vector<int> &end);
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
  float2DReg(const std::shared_ptr<float2DReg> old, const int iax1,
             const bool rev1, const int iax2, const bool rev2,
             const std::vector<int> &ipos, const std::vector<int> &beg,
             const std::vector<int> &end);
  /*!
    Allocate data for vector
    */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new float2D(boost::extents[ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
    Create a window of the current vector
    \param nw,fw,jw Window parameters
    */
  std::shared_ptr<float2DReg> window(const std::vector<int> &nw,
                                     const std::vector<int> &fw,
                                     const std::vector<int> &jw) const;
  /*!
   Make a copy of the vector
   */
  std::shared_ptr<float2DReg> clone() const;
  /*!
 Make a copy of the vector space
 */
  std::shared_ptr<float2DReg> cloneSpace() const;
  /*!
 Deallocate storage for vector, turn into vector space
  */
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<float2D> _mat;  ///< Storage for vector

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
  void initData(std::shared_ptr<SEP::hypercube> hyper, const float2D &vals);
};
}  // namespace SEP
#endif
