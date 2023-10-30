#ifndef JSON_PARAM_OBJ
#define JSON_PARAM_OBJ 1
#include <memory>
#include "json.h"
#include "paramObj.h"
namespace SEP {
/*!
  Get arguments from JSON parameters
  */
class jsonParamObj : public paramObj {
 public:
  /*!
    Initalize JSON object from JSON  value
    \param inp JSON value
    */
  jsonParamObj(Json::Value &inp);
  /*!
Get an integer from a file

\param arg Name of the prameter
*/
  virtual int getInt(const std::string &arg) const;
  /*!
Get an integer from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual int getInt(const std::string &arg, const int def) const;
  /*!
  Get a float from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
 */
  virtual float getFloat(const std::string &, const float def) const;
  /*!
Get a float from a file

\param arg Name of the prameter
*/
  virtual float getFloat(const std::string &) const;
  /*!
Get a string  from a file

\param arg Name of the prameter
*/
  virtual std::string getString(const std::string &arg) const;
  /*!
Get a string from a file

\param tag Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::string getString(const std::string &arg,
                                const std::string &def) const;
  /*!
Get boolean from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual bool getBool(const std::string &, const bool def) const;
  /*!
Get a boolean from a file

\param arg Name of the prameter
*/
  virtual bool getBool(const std::string &) const;
  /*!
Get integer from a file

\param arg Name of the prameter
\param nvals Number of values
*/
  virtual std::vector<int> getInts(const std::string &arg,
                                   const int nvals) const;
  /*!
Get integer from a file

\param arg Name of the prameter
\param defs Default values
*/
  virtual std::vector<int> getInts(const std::string &arg,
                                   const std::vector<int> &defs) const;
  /*!
  Get an floats from a file

  \param arg Name of the prameter
  \param nval Number of values to look for
 */
  virtual std::vector<float> getFloats(const std::string &arg,
                                       const int nvals) const;
  /*!
Get floats from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::vector<float> getFloats(const std::string &arg,
                                       const std::vector<float> &defs) const;
  /*!
 Output a message
 \param err Message to output
 */
  virtual void message(const std::string &err) const;
  /*!
Output a message and exit with an error
\param err Message to output
*/
  virtual void error(const std::string &err) const;
  /*!
   Add additional parameters

   \param pars List of parameters to add
   */
  virtual void addParams(std::map<std::string, std::string> &pars);

 protected:
  Json::Value jsonArgs;  ///< JSON parameter object
};

}  // namespace SEP

#endif