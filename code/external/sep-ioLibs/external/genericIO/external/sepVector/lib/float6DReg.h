#ifndef float6d_reg_h
#define float6d_reg_h 1
#include <floatHyper.h>
#include "boost/multi_array.hpp"


typedef boost::multi_array<float, 6> float6D;
namespace SEP {
/*!
A regular sampled 6-D function with float storage
*/
class float6DReg : public floatHyper {
 public:
  /*!
   Create a 6-D float vector from a hypercube
        \param hyper Hypercube describing RSF

   */
  float6DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
 Create a 6-D float vector from just lengths
      \param n1,n2,n3,n4,n5,n6 Dimensions of the hypercube

 */
  float6DReg(const int n1, const int n2, const int n3, const int n4,
             const int n5, const int n6) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    a.push_back(SEP::axis(n4));
    a.push_back(SEP::axis(n5));
    a.push_back(SEP::axis(n6));

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
Create a 6-D float vector from axes
    \param a1,a2,a3,a4,a5,a6 Axes if the hypercube

*/
  float6DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
             const SEP::axis &a4, const SEP::axis &a5, const SEP::axis &a6) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    a.push_back(a4);
    a.push_back(a5);
    a.push_back(a6);

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 6-D float vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Values to fill vector with

 */
  float6DReg(std::shared_ptr<SEP::hypercube> hyper, const float6D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 6-D float vector from just lengths
    \param n1,n2,n3,n4,n5,n6 Dimensions of the hypercube
    \param vals Values to fill vector with

*/
  float6DReg(const int n1, const int n2, const int n3, const int n4,
             const int n5, const int n6, float6D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    a.push_back(SEP::axis(n4));
    a.push_back(SEP::axis(n5));
    a.push_back(SEP::axis(n6));

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
 Create a 6-D float vector from axes
      \param a1,a2,a3,a4,a5,a6 Axes if the hypercube
      \param vals Values to fill vector with
 */
  float6DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
             const SEP::axis &a4, const SEP::axis &a5, const SEP::axis &a6,
             float6D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    a.push_back(a4);
    a.push_back(a5);
    a.push_back(a6);

    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
     Allocate data for vector
     */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(
        new float6D(boost::extents[ns[5]][ns[4]][ns[3]][ns[2]][ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
Make a copy of the vector
*/
  std::shared_ptr<float6DReg> clone() const;
  /*!
 Make a copy of the vector space
 */
  std::shared_ptr<float6DReg> cloneSpace() const;
  /*!
                                  Return a subset of the vector
                                  \param nw,fw,jw Windowing parameters
                                  */
  std::shared_ptr<float6DReg> window(const std::vector<int> &nw,
                                     const std::vector<int> &fw,
                                     const std::vector<int> &jw) const;
  /*!
  Deallocate storage for vector, turn into vector space
   */
  virtual void cleanMemory() {
    _mat = 0;
    setSpace();
  }
  std::shared_ptr<float6D> _mat;  ///< Storage for vector

 protected:
  /*!
  Initialize vector without providing data (zero)
  \param hyp Hypercube describing dataset
*/
  void initNoData(std::shared_ptr<SEP::hypercube> hyp);
  /*!
Initialize with data
\param hyper Hypercube describing space
\param vals Data to copy in
*/
  void initData(std::shared_ptr<SEP::hypercube> hyper, const float6D &vals);
};
}  // namespace SEP
#endif
