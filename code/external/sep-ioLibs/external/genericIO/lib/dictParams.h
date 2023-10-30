#ifndef DICT_PARAM_FUNC_H
#define DICT_PARAM_FUNC_H 1
#include <stdbool.h>
#include <map>
#include <string>
#include "paramObj.h"
namespace SEP {
/*!
  Python parameter object.
    Emulate command line parameter reading from a python map
  */
class dictParams : public SEP::paramObj {
 public:
  /*!
     Initialize python parameter object
     \param pars Parameters in a dictionary -> map
     */
  dictParams(std::map<std::string, std::string> pars);

  /*!
  Reset the parameter database

\param pars Name of the prameter
*/
  void resetParams(const std::map<std::string, std::string> pars);

  /*!
  Add to parameter database

\param pars Name of the prameter
*/
  void addParams(const std::map<std::string, std::string> pars);

  /*!
Get an integer from a file

\param arg Name of the prameter
*/
  virtual int getInt(const std::string &arg) const override;
  /*!
Get an integer from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual int getInt(const std::string &arg, const int def) const override;
  /*!
  Get a float from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
 */
  virtual float getFloat(const std::string &, const float def) const override;
  /*!
Get a float from a file

\param arg Name of the prameter
*/
  virtual float getFloat(const std::string &) const override;
  /*!
Get a string  from a file

\param arg Name of the prameter
*/
  virtual std::string getString(const std::string &arg) const override;
  /*!
Get a string from a file

\param tag Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::string getString(const std::string &arg,
                                const std::string &def) const override;
  /*!
Get boolean from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual bool getBool(const std::string &, const bool def) const override;
  /*!
Get a boolean from a file

\param arg Name of the prameter
*/
  virtual bool getBool(const std::string &) const override;
  /*!
Get integer from a file

\param arg Name of the prameter
\param nvals Number of values
*/
  virtual std::vector<int> getInts(const std::string &arg,
                                   const int nvals) const override;
  /*!
Get integer from a file

\param arg Name of the prameter
\param defs Default values
*/
  virtual std::vector<int> getInts(const std::string &arg,
                                   const std::vector<int> &defs) const override;
  /*!
  Get an floats from a file

  \param arg Name of the prameter
  \param nval Number of values to look for
 */
  virtual std::vector<float> getFloats(const std::string &arg,
                                       const int nvals) const override;
  /*!
Get floats from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::vector<float> getFloats(
      const std::string &arg, const std::vector<float> &defs) const override;
  /*!
Output a message and exit with an error
\param err Message to output
*/
  virtual void error(const std::string &errm) const override;
  /*!
 Output a message
 \param err Message to output
 */
  virtual void message(const std::string &msg) const override;
  /*!
    Split a string into parts
    \param in String to split
    */
  std::vector<std::string> splitString(const std::string &in) const;

  /*!
   Add additional parameters

   \param pars List of parameters to add
   */
  virtual void addParams(std::map<std::string, std::string> &pars)override;

 private:
  std::map<std::string, std::string>
      _pars;  ///< Map containing the equivielent of command line arguments
};
}  // namespace SEP

#endif
