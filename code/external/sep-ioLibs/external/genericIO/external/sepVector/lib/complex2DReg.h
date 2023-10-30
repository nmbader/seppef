#ifndef complex2d_reg_h
#define complex2d_reg_h 1
#include <complexHyper.h>
#include "boost/multi_array.hpp"
typedef boost::multi_array<std::complex<float>, 2> complex2D;
namespace SEP {
/*!
A regular sampled 2-D function with complex float storage
*/
class complex2DReg : public complexHyper {
 public:
  /*!
   Create a 2-D complex float vector from a hypercube
        \param hyper Hypercube describing RSF

   */
  complex2DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 2-D complex float vector from just lengths
    \param n1,n2 Dimensions of the hypercube

*/
  complex2DReg(const int n1, const int n2) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 2-D complex float vector from axes
      \param a1,a2 Axes if the hypercube

 */
  complex2DReg(const SEP::axis &a1, const SEP::axis &a2) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 2-D complex float vector from a hypercube
      \param hyper Hypercube describing RSF
      \param vals Data to fill vector with

 */
  complex2DReg(std::shared_ptr<SEP::hypercube> hyper, const complex2D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 2-D complex float vector from just lengths
  \param n1,n2 Dimensions of the hypercube
  \param vals Data to fill vector with

*/
  complex2DReg(int n1, int n2, complex2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(SEP::axis(n1));
    a.push_back(SEP::axis(n2));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
 Create a 2-D complex float vector from axes
      \param a1,a2 Axes if the hypercube
      \param vals Data to fill vector with
 */
  complex2DReg(SEP::axis &a1, SEP::axis &a2, complex2D &vals) {
    std::vector<SEP::axis> a;
    a.push_back(a1);
    a.push_back(a2);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
    Allocate data for vector
    */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new complex2D(boost::extents[ns[1]][ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
   Make a copy of the vector
   */
  std::shared_ptr<complex2DReg> clone() const;
  /*!
 Make a copy of the vector space
 */
  std::shared_ptr<complex2DReg> cloneSpace() const;
  /*!
   Create a window of the current vector
   \param nw,fw,jw Window parameters
   */
  std::shared_ptr<complex2DReg> window(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw) const;
  /*!
 Deallocate storage for vector, turn into vector space
  */
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<complex2D> _mat;  ///< Storage for vector

 private:
  /*!
Initialize without data
\param hyper Hypercube describing space
*/
  void initNoData(std::shared_ptr<SEP::hypercube> hyp);
  /*!
    Initialize with data
    \param hyper Hypercube describing space
    \param vals Data to copy in
 */
  void initData(std::shared_ptr<SEP::hypercube> hyp, const complex2D &vals);
};
}  // namespace SEP
#endif
