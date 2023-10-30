#ifndef RECTFILTER_H
#define RECTFILTER_H 1
#include "float1DReg.h"
#include "float2DReg.h"
#include "float3DReg.h"
#include "float4DReg.h"
#include "float5DReg.h"
namespace SEP {
/*!
  A rectangular 1-D filter

   These filter describe everything as n-dimensional box. This is done
   for code efficiency to avoid indirection with lags. It handles irregular
   filters by zeroing updates to filter coeficients that aren't allowed to
  change.

*/
class rectFilter1D : public SEP::float1DReg {
 public:
  /*!
  Create a rectangular filter from a box
    \param box A vector describing the rectangular space of the filter
    \param beg Location of the fixed (1) coefficient [only valid pef]
    \param pef Whether or not this a Prediction Error Filter
*/
  rectFilter1D(const std::vector<int> &box, const std::vector<int> &beg,
               bool pef = false);
  /*!
 Clone the filter
*/
  std::shared_ptr<rectFilter1D> clone() const;
  /*!
Clone the filter vector space
*/
  std::shared_ptr<rectFilter1D> cloneSpace() const;
  /*!
Ignore updates to the rectangular filter that aren't part of the irregular
filter
*/
  void cleanFilter() {
    zeroNonCoefs();
    (*_mat)[_f[0]] = 1;
  }
  /*!
     Zdero the coeficients that aren't allowed to chabnge
     */
  void zeroNonCoefs();
  /*!
     Return whether or not the filter is a PEF
     */
  bool isPef() { return _pef; }
  /*!
    Set filter basics

    \param box - Description of box containing filter
    \param f   - Location of 0 lag coefficient
    \param pef - Whether or not a boolean
  */

  void setBasics(const std::vector<int> &box, const std::vector<int> &f,
                 const bool pef) {
    _n = box;
    _f = f;
    for (auto i = 0; i < box.size(); i++) {
      _e.push_back(_n[i] - _f[i] - 1);
    }
    _pef = pef;
  }

 public:
  std::vector<int> _n;  ///< Rectangular filter size
  std::vector<int> _f;  ///< Offset location of 1 from top-left
  std::vector<int> _e;  ///< Offset location of 1 from bottom-right
  bool _pef;            ///< Whether or not a boolean
};
/*!
  A rectangular 2-D filter

     These filter describe everything as n-dimensional box. This is done
   for code efficiency to avoid indirection with lags. It handles irregular
   filters by zeroing updates to filter coeficients that aren't allowed to
  change.
*/
class rectFilter2D : public SEP::float2DReg {
 public:
  /*!
Create a rectangular filter from a box
  \param box A vector describing the rectangular space of the filter
  \param beg Location of the fixed (1) coefficient [only valid pef]
  \param pef Whether or not this a Prediction Error Filter
*/
  rectFilter2D(const std::vector<int> &box, const std::vector<int> &beg,
               bool pef = false);
  /*!
Clone the filter
*/
  std::shared_ptr<rectFilter2D> clone() const;
  /*!
Clone the filter vector space
*/
  std::shared_ptr<rectFilter2D> cloneSpace() const;

  /*!
    Zdero the coeficients that aren't allowed to chabnge
    */
  void zeroNonCoefs();

  /*!
Ignore updates to the rectangular filter that aren't part of the irregular
filter
*/
  void cleanFilter() {
    zeroNonCoefs();
    (*_mat)[_f[1]][_f[0]] = 1;
  }
  /*!
   Return whether or not the filter is a PEF
   */
  bool isPef() { return _pef; }
  /*!
   Set the basic properites of the filter
   \param box Describe the rectangular nature of the filter
   \param f Location of the 1 element
   \param pef Whther or not a PEF
   */
  void setBasics(const std::vector<int> &box, const std::vector<int> &f,
                 const bool pef) {
    _n = box;
    _f = f;
    for (auto i = 0; i < box.size(); i++) {
      _e.push_back(_n[i] - _f[i] - 1);
    }
    _pef = pef;
  }

 public:
  std::vector<int> _n;  ///< Rectangular filter size
  std::vector<int> _f;  ///< Offset location of 1 from top-left
  std::vector<int> _e;  ///< Offset location of 1 from bottom-right
  bool _pef;            ///< Whether or not a boolean
};
/*!
  A rectangular 3-D filter

     These filter describe everything as n-dimensional box. This is done
   for code efficiency to avoid indirection with lags. It handles irregular
   filters by zeroing updates to filter coeficients that aren't allowed to
  change.
*/
class rectFilter3D : public SEP::float3DReg {
 public:
  /*!
Create a rectangular filter from a box
  \param box A vector describing the rectangular space of the filter
  \param beg Location of the fixed (1) coefficient [only valid pef]
  \param pef Whether or not this a Prediction Error Filter
*/
  rectFilter3D(const std::vector<int> &box, const std::vector<int> &beg,
               bool pef = false);
  /*!
Clone the filter
*/
  std::shared_ptr<rectFilter3D> clone() const;
  /*!
Clone the filter vector space
*/
  std::shared_ptr<rectFilter3D> cloneSpace() const;
  /*!
    Zdero the coeficients that aren't allowed to chabnge
    */
  void zeroNonCoefs();
  /*!
Ignore updates to the rectangular filter that aren't part of the irregular
filter
*/
  void cleanFilter() {
    zeroNonCoefs();
    (*_mat)[_f[2]][_f[1]][_f[0]] = 1;
  }
  /*!
   Return whether or not the filter is a PEF
   */
  bool isPef() { return _pef; }
  /*!
 Set the basic properites of the filter
 \param box Describe the rectangular nature of the filter
 \param f Location of the 1 element
 \param pef Whther or not a PEF
 */
  void setBasics(const std::vector<int> &box, const std::vector<int> &f,
                 const bool pef) {
    _n = box;
    _f = f;
    for (auto i = 0; i < box.size(); i++) {
      _e.push_back(_n[i] - _f[i] - 1);
    }
    _pef = pef;
  }

 public:
  std::vector<int> _n;  ///< Rectangular filter size
  std::vector<int> _f;  ///< Offset location of 1 from top-left
  std::vector<int> _e;  ///< Offset location of 1 from bottom-right
  bool _pef;            ///< Whether or not a boolean
};
/*!
  A rectangular 4-D filter

     These filter describe everything as n-dimensional box. This is done
   for code efficiency to avoid indirection with lags. It handles irregular
   filters by zeroing updates to filter coeficients that aren't allowed to
  change.
*/
class rectFilter4D : public SEP::float4DReg {
 public:
  /*!
Create a rectangular filter from a box
  \param box A vector describing the rectangular space of the filter
  \param beg Location of the fixed (1) coefficient [only valid pef]
  \param pef Whether or not this a Prediction Error Filter
*/
  rectFilter4D(const std::vector<int> &box, const std::vector<int> &beg,
               bool pef = false);
  /*!
Clone the filter
*/
  std::shared_ptr<rectFilter4D> clone() const;
  /*!
Clone the filter vector space
*/
  std::shared_ptr<rectFilter4D> cloneSpace() const;
  /*!
    Zdero the coeficients that aren't allowed to chabnge
    */
  void zeroNonCoefs();
  /*!
Ignore updates to the rectangular filter that aren't part of the irregular
filter
*/
  void cleanFilter() {
    zeroNonCoefs();
    (*_mat)[_f[3]][_f[2]][_f[1]][_f[0]] = 1;
  }
  /*!
 Return whether or not the filter is a PEF
 */
  bool isPef() { return _pef; }
  /*!
 Set the basic properites of the filter
 \param box Describe the rectangular nature of the filter
 \param f Location of the 1 element
 \param pef Whther or not a PEF
 */
  void setBasics(const std::vector<int> &box, const std::vector<int> &f,
                 const bool pef) {
    _n = box;
    _f = f;
    for (auto i = 0; i < box.size(); i++) {
      _e.push_back(_n[i] - _f[i] - 1);
    }
    _pef = pef;
  }

 public:
  std::vector<int> _n;  ///< Rectangular filter size
  std::vector<int> _f;  ///< Offset location of 1 from top-left
  std::vector<int> _e;  ///< Offset location of 1 from bottom-right
  bool _pef;            ///< Whether or not a boolean
};

}  // namespace SEP
#endif