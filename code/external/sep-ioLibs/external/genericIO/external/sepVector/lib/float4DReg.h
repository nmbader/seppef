#ifndef float4d_reg_h
#define float4d_reg_h 1
#include <floatHyper.h>
#include "boost/multi_array.hpp"


typedef boost::multi_array<float, 4> float4D;
namespace SEP {
/*!
A regular sampled 4-D function with float storage
*/
class float4DReg : public floatHyper {
 public:
  /*!
  Default object
  */
  float4DReg() { ; }
  /*!
   Create a 4-D float vector from a hypercube
   \param hyper Hypercube describing RSF
   */
  float4DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 4-D float vector from just lengths
  \param n1,n2,n3,n4 Dimensions of the hypercube

*/
  float4DReg(const int n1, const int n2, const int n3, const int n4) {
	  std::cerr<<"what is going on 1"<<std::endl;
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    a.push_back(SEP::axis(n4));
	  std::cerr<<"2hat is going on 1"<<std::endl;
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
	  std::cerr<<"3hat is going on 1"<<std::endl;
    initNoData(hyp);
	  std::cerr<<"4hat is going on 1"<<std::endl;
  }
  /*!
   Create a 4-D float vector from axes
        \param a1,a2,a3,a4 Axes if the hypercube

   */
  float4DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
             const SEP::axis &a4) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    a.push_back(a4);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
Create a 4-D float vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Vaules for the vector

*/
  float4DReg(std::shared_ptr<SEP::hypercube> hyper, const float4D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 4-D float vector from a hypercube
      \param n1,n2,n3,n4 length of axes for hypercube
      \param vals Vaules for the vector

*/
  float4DReg(const int n1, const int n2, const int n3, const int n4,
             float4D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    a.push_back(SEP::axis(n4));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
 Create a -D float vector from axes
      \param a1,a2,a3,a4 Axes if the hypercube
      \param vals Vaules for the vector
 */
  float4DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
             const SEP::axis &a4, float4D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    a.push_back(a4);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
      Allocate data for vector
      */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new float4D(boost::extents[ns[3]][ns[2]][ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
Make a copy of the vector
*/
  std::shared_ptr<float4DReg> clone() const;
  /*!
  Make a copy of the vector space
  */
  std::shared_ptr<float4DReg> cloneSpace() const;
  /*!
    Return a subset of the vector
    \param nw,fw,jw Windowing parameters
    */
  std::shared_ptr<float4DReg> window(const std::vector<int> &nw,
                                     const std::vector<int> &fw,
                                     const std::vector<int> &jw) const;
  /*!
  Deallocate storage for vector, turn into vector space
   */
  virtual void cleanMemory() {
    _mat = 0;
    setSpace();
  }
  std::shared_ptr<float4D> _mat;  ///< Storage for vector

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
  void initData(std::shared_ptr<SEP::hypercube> hyper, const float4D &vals);
};
}  // namespace SEP
#endif
