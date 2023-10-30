#ifndef HYPERCUBE_H
#define HYPERCUBE_H 1
#include <axis.h>
#include <memory>
#include <sstream>
#include <vector>
namespace SEP {
/*!
Describing a multi-dimensional regular sampled function
*/
class hypercube {
 public:
  hypercube() {}
  //! Create a hypercube from another hyperxube
  /*!
    \param hyper Create a hypercube from another hypercube
  */
  hypercube(const hypercube &hyper);
  /// Override equal
  hypercube &operator=(const hypercube &t) {
    axes = t.getAxes();
    n123 = t.getN123();
    return *this;
  }
  //! Create a hypercube of ndim dimensions
  /*!
    \param ndim Number of diemensions
  */
  hypercube(const int ndim);
  //! Create a hypercube from a shared pointer hypercube
  /*!
    \param hyper Shared pointer hypercube
  */
  hypercube(const std::shared_ptr<hypercube> hyper);

  //! Create a hypercube from a single axis
  /*!
    \param a1 First axis
  */
  hypercube(const SEP::axis &a1) {
    std::vector<SEP::axis> as;
    as.push_back(a1);
    setAxes(as);
  }

  //! Create a hypercube from two axes
  /*!
    \param a1 First axis
    \param a2 Second axis
  */
  hypercube(const SEP::axis &a1, const SEP::axis &a2) {
    std::vector<SEP::axis> as;
    as.push_back(a1);
    as.push_back(a2);
    setAxes(as);
  }
  //! Create a hypercube from three axes
  /*!
    \param a1 First axis
    \param a2 Second axis
    \param a3 Third axis
  */
  hypercube(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3) {
    std::vector<SEP::axis> as;
    as.push_back(a1);
    as.push_back(a2);
    as.push_back(a3);
    setAxes(as);
  }

  //! Create a hypercube from three axes
  /*!
    \param axes A list of axes to create hypercube from
  */
  hypercube(const std::vector<SEP::axis> &axes);

  //! Clone a hypercube from another hypercube

  std::shared_ptr<hypercube> clone() const {
    std::shared_ptr<hypercube> a(new hypercube(*this));
    return a;
  }
  //! Create a hypercube from four axes
  /*!
    \param a1 First axis
    \param a2 Second axis
    \param a3 Third axis
    \param a4 Fourth axis
  */
  hypercube(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
            const SEP::axis &a4) {
    std::vector<SEP::axis> as;
    as.push_back(a1);
    as.push_back(a2);
    as.push_back(a3);
    as.push_back(a4);
    setAxes(as);
  }
  //! Create a hypercube from five axes
  /*!
    \param a1 First axis
    \param a2 Second axis
    \param a3 Third axis
    \param a4 Fourth axis
    \param a5 Fifth axis
  */
  hypercube(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
            const SEP::axis &a4,const SEP::axis &a5) {
    std::vector<SEP::axis> as;
    as.push_back(a1);
    as.push_back(a2);
    as.push_back(a3);
    as.push_back(a4);
    as.push_back(a5);
    setAxes(as);
  }
  //! Create a hypercube from six axes
  /*!
    \param a1 First axis
    \param a2 Second axis
    \param a3 Third axis
    \param a4 Fourth axis
    \param a5 Fifth axis
    \param a6 Sixth axis
  */
  hypercube(const SEP::axis &a1, const SEP::axis &a2, const SEP::axis &a3,
            const SEP::axis &a4,const SEP::axis &a5,const SEP::axis &a6) {
    std::vector<SEP::axis> as;
    as.push_back(a1);
    as.push_back(a2);
    as.push_back(a3);
    as.push_back(a4);
    as.push_back(a5);
    as.push_back(a6);
    setAxes(as);
  }

  //! Check to see if two hypercubes are the same
  /*!
    \param hyper Hypercube to compare to
  */
  bool checkSame(const std::shared_ptr<hypercube> hyper) const;
  //! Set axes
  /*!
    \param axes List of axis
  */
  void setAxes(const std::vector<SEP::axis> &axes);
  //! Set axis
  /*!
    \param idim - Axis to set (starting with 1)
    \param ax  - Axis
  */
  void setAxis(const int idim, const SEP::axis &ax);
  //! Get an axis
  /*!
    \param idim Return a given axis
  */
  SEP::axis getAxis(const int idim) const;
  //! Get total number of elements

  long long getN123() const { return n123; }
  //! Add to a stream information about hypercube
  /*!
    \param stream Stream to add to
  */
  void infoStream(std::stringstream &stream);
  //! Get the ns from the hypercube

  std::vector<int> getNs() const;

  //! Shrink the number of dimensions of the hypercube
  /*!
    \param nmax to shrink the hypercube to
  */
  void shrinkDimension(const int nmax) { axes.resize(nmax); }

  //! Append axis
  /*!
    \param ax to append
  */
  void addAxis(axis &ax) { axes.push_back(ax); }
  //!  Remove axes
  void deallocate() { axes.clear(); }

  //! Delete hypercube object
  ~hypercube() { this->deallocate(); }

  //! Return up to nmax axes from the hypercube
  /*!
    \param nmax  Maximum number of axes to return
  */
  std::vector<SEP::axis> returnAxes(const int nmax) const;
  //! Return the number of dimensions in the hypercube
  int getNdim() const { return (int)axes.size(); }
  //! Return the number of dimensions in the hypercube greater than 1 in size
  int getNdimG1() const;
  //! Return all the axes in hypercube
  std::vector<SEP::axis> getAxes() const;
  //! Return all the axes in hypercube up to nmin, pad if 1
  /*!
  \param nmin  Minimum number of axes to return
*/
  std::vector<SEP::axis> getAxes(const int nmin) const;
  //! Check to see if hypercube is the same size
  /*!
  \param other Hypercube to compare to
*/
  bool sameSize(const hypercube &other) const;
  //! Check to see if hypercube is the same size
  /*!
  \param other Hypercube to compare to
*/
  bool sameSize(const std::shared_ptr<hypercube> &other) const;

 private:
  //! Internal function to initialize a hypercube object
  /*!
    \param axes  Axes to initialize with
  */
  void initNd(const std::vector<SEP::axis> &axes);

 protected:
  long long n123;               ///< Number of samples
  std::vector<SEP::axis> axes;  ///< List of axes
};
}  // namespace SEP
#endif
