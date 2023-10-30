#ifndef AXIS_H
#define AXIS_H 1
#include <sstream>
#include <string>

namespace SEP {

/*!
  A class describing  a regular sampled axis.
*/
class axis {
 public:
  axis() {}  /// Default constructor does nothing

  //! Create an axis object
  /*!
    \param n Number of samples in the axis
    \param o Origin of axis
    \param d Sampling between samples
    \param label Label for the axis
    \param unit Unit for axis
  */
  axis(const int n, float o = 0., float d = 1.,
       std::string label = std::string(), std::string unit = std::string());

  //! Assign axis from another axis
  axis &operator=(const axis &t) {
    n = t.n;
    o = t.o;
    d = t.d;
    label = t.label;
    unit = t.unit;
    return *this;
  }
  //! Add to info stream
  /*!
    \param stream Stream to add info about axis to
  */
  void infoStream(std::stringstream &stream);

  //! Is the passed along axis represent the same axis
  /*!
    \param ax Axis to compare object to
  */
  bool same_axis(const axis &ax) const;
  int n;              ///< Number of samples in axis
  float o;            ///< Origin os axis
  float d;            ///< Sampling of axis
  std::string label;  ///< Label for axis
  std::string unit;   ///< Unit for axis

  ~axis() {}

 private:
};  // namespace SEP
}  // namespace SEP
#endif
