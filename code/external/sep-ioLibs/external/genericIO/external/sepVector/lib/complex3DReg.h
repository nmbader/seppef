#ifndef complex3d_reg_h
#define complex3d_reg_h 1
#include <complexHyper.h>
#include "boost/multi_array.hpp"
typedef boost::multi_array<std::complex<float>, 3> complex3D;
//! Namespace for SEP code
namespace SEP {
/*!
A regular sampled 3-D function with complex float storage
*/
class complex3DReg : public complexHyper {
 public:
  /*!
   Create a 3-D complex float vector from a hypercube
        \param hyper Hypercube describing RSF

   */
  complex3DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 3-D complex float vector from just lengths
     \param n1,n2,n3 Dimensions of the hypercube

*/
  complex3DReg(const int n1, const int n2, const int n3) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 3-D complex float vector from axes
      \param a1,a2,a3 Axes if the hypercube

 */
  complex3DReg(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 3-D complex float vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Vaule to file vector with
 */
  complex3DReg(std::shared_ptr<SEP::hypercube> hyper, const complex3D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 3-D complex float vector from just lengths
  \param n1,n2,n3 Dimensions of the hypercube
  \param vals Vaule to file vector with

*/
  complex3DReg(int n1, int n2, int n3, complex3D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    a.push_back(SEP::axis(n3));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
Create a 3-D complex float vector from axes
    \param a1,a2,a3 Axes if the hypercube
    \param vals Vaule to file vector with

*/
  complex3DReg(SEP::axis &a1, SEP::axis &a2, SEP::axis &a3, complex3D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    a.push_back(a3);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
  Allocate data for vector
  */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new complex3D(boost::extents[ns[2]][ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
   Make a copy of the vector
   */
  std::shared_ptr<complex3DReg> clone() const;
  /*!
Make a copy of the vector space
*/
  std::shared_ptr<complex3DReg> cloneSpace() const;
  /*!
  Deallocate storage for vector, turn into vector space
   */
  virtual void cleanMemory() {
    _mat = 0;
    setSpace();
  }
  /*!
    Return a subset of the vector
    \param nw,fw,jw Windowing parameters
    */
  std::shared_ptr<complex3DReg> window(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw) const;
  std::shared_ptr<complex3D> _mat;  ///< Storage for vector

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

  void initData(std::shared_ptr<SEP::hypercube> hyper, const complex3D &vals);
  /*!
Initialize with data
\param hyper Hypercube describing space
\param vals Data to copy in
*/

  void initData(std::shared_ptr<SEP::hypercube> hyper,
                std::shared_ptr<complex3D> vals);
};
}  // namespace SEP
#endif
