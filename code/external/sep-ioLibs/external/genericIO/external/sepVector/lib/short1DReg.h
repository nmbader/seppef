#ifndef short1d_reg_h
#define short1d_reg_h 1
#include <shortHyper.h>
#include <cstdint>
#include <iostream>
#include "boost/multi_array.hpp"
typedef boost::multi_array<short, 1> short1D;
namespace SEP {
/*!
A regular sampled 1-D function with short storage
*/
class short1DReg : public shortHyper {
 public:
  /*!
 Create a 1-D short vector from a hypercube
      \param hyper Hypercube describing RSF

 */
  short1DReg(std::shared_ptr<SEP::hypercube> hyper) { initNoData(hyper); }
  /*!
Create a 1-D short vector from just lengths
    \param n  Dimensions of the hypercube

*/
  short1DReg(const int n) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initNoData(hyp);
  }
  /*!
 Create a 1-D short vector from axes
      \param a1 Axes of the hypercube

 */
  short1DReg(const SEP::axis &a1) {
    std::vector<SEP::axis> as(1, a1);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    initNoData(hyp);
  }
  /*!
 Create a 1-D short vector from just lengths
      \param hyper Hypercube describing RSF
     \param vals Values to fill vector with

 */
  short1DReg(std::shared_ptr<SEP::hypercube> hyper, const short1D &vals) {
    initData(hyper, vals);
  }
  /*!
Create a 1-D short vector from just lengths
  \param n  Dimensions of the hypercube
 \param vals Values to fill vector with

*/
  short1DReg(const int n, short1D &vals) {
    std::vector<SEP::axis> a(1, SEP::axis(n));
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(a));
    initData(hyp, vals);
  }
  /*!
Create a 1-D short vector from axes
     \param ax Axes of the hypercube
    \param vals Values to fill vector with

*/
  short1DReg(const SEP::axis &ax, const short1D &vals) {
    std::vector<SEP::axis> as(1, ax);
    std::shared_ptr<SEP::hypercube> hyp(new SEP::hypercube(as));
    initData(hyp, vals);
  }
  /*!
   Allocate data for vector
   */
  void allocate() {
    std::vector<int> ns = getHyper()->getNs();
    _mat.reset(new short1D(boost::extents[ns[0]]));

    setData(_mat->data());
    ;
  }
  /*!
 Make a copy of the vector
 */
  std::shared_ptr<short1DReg> clone() const;
  /*!
Make a copy of the vector
*/
  std::shared_ptr<short1DReg> cloneSpace() const;
  /*!
Deallocate storage for vector, turn into vector space
 */
  virtual void cleanMemory() { setSpace(); }
  std::shared_ptr<short1D> _mat;  ///< Storage for vector

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
  void initData(std::shared_ptr<SEP::hypercube> hyper, const short1D &vals);
};

}  // namespace SEP
#endif
