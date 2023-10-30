#ifndef PARAM_OBJ_H
#define PARAM_OBJ_H 1
#include <stdbool.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "SEPException.h"

namespace SEP {
/*!
  Basic parameter object (virtual)
  */
class paramObj {
 public:
  std::string _type;  ///< Parameter type

  /*!
     Default param object
     */
  paramObj() { _type = "NONE"; }
  /*!
Output a message and exit with an error
\param err Message to output
*/
  virtual void error(const std::string& err) const { throw SEPException(err); }
  /*!
   Output a message
   \param err Message to output
   */

  virtual void message(const std::string& msg) const {
    std::cerr << msg << std::endl;
  }
  /*!
Get an integer from a file

\param arg Name of the prameter
*/
  virtual int getInt(const std::string& arg) const {
    if (arg == "")
      ;
    return 0;
  }
  /*!
Get an integer from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual int getInt(const std::string& arg, const int def) const {
    if (arg == "" || def == 0)
      ;
    return 0;
  }
  /*!
  Get a float from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
 */
  virtual float getFloat(const std::string& arg, const float def) const {
    if (arg == "" || def == 0)
      ;
    return 0.;
  }
  /*!
Get a float from a file

\param arg Name of the prameter
*/
  virtual float getFloat(const std::string& arg) const {
    if (arg == "")
      ;
    return 1.;
  }
  /*!
Get a string  from a file

\param arg Name of the prameter
*/
  virtual std::string getString(const std::string& arg) const {
    if (arg == "")
      ;
    return std::string("");
  }
  /*!
Get a string from a file

\param tag Name of the prameter
\param def Default value (if not found in file)
*/
  virtual std::string getString(const std::string& arg,
                                const std::string& def) const {
    if (arg == def)
      ;
    return std::string("");
  }
  /*!
Get boolean from a file

\param arg Name of the prameter
\param def Default value (if not found in file)
*/
  virtual bool getBool(const std::string& arg, const bool def) const {
    if (arg == "" || def == 0)
      ;
    return false;
  }
  /*!
  Get a boolean from a file

  \param arg Name of the prameter
  */
  virtual bool getBool(const std::string& arg) const {
    if (arg == "")
      ;
    return false;
  }
  /*!
  Get integer from a file

  \param arg Name of the prameter
  \param nvals Number of values
  */
  virtual std::vector<int> getInts(const std::string& arg,
                                   const int nvals) const {
    if (arg == "" || nvals == 0)
      ;
    return std::vector<int>(1);
  }
  /*!
  Get integer from a file

  \param arg Name of the prameter
  \param defs Default values
  */
  virtual std::vector<int> getInts(const std::string& arg,
                                   const std::vector<int>& defs) const {
    if (arg == "" && defs.size() == -1)
      ;
    return std::vector<int>(1);
  }
  /*!
  Get an floats from a file

  \param arg Name of the prameter
  \param nval Number of values to look for
  */
  virtual std::vector<float> getFloats(const std::string& arg,
                                       int nvals) const {
    if (arg == "" || nvals == 0)
      ;
    return std::vector<float>(1);
  }
  /*!
  Get floats from a file

  \param arg Name of the prameter
  \param def Default value (if not found in file)
  */
  virtual std::vector<float> getFloats(const std::string& arg,
                                       const std::vector<float>& defs) const {
    if (arg == "" || defs.size() == -1)
      ;
    return std::vector<float>(1);
  }

  /*!
  Get an integer from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter
  */
  virtual int getDocInt(const std::string& arg, const std::string& doc);
  /*!
  Get an integer from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param def Default value (if not found in file)
  */
  virtual int getDocInt(const std::string& arg, const std::string& doc,
                        const int def);
  /*!
  Get a float from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param def Default value (if not found in file)
  */
  virtual float getDocFloat(const std::string& arg, const std::string& doc,
                            const float def);
  /*!
  Get a float from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  */
  virtual float getDocFloat(const std::string& arg, const std::string& doc);
  /*!
  Get a string  from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  */
  virtual std::string getDocString(const std::string& arg,
                                   const std::string& doc);
  /*!
  Get a string from a file

  \param tag Name of the prameter

  \param doc Documentation for parameter

  \param def Default value (if not found in file)
  */
  virtual std::string getDocString(const std::string& arg,
                                   const std::string& doc,
                                   const std::string& def);

  /*!
  Get boolean from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param def Default value (if not found in file)
  */
  virtual bool getDocBool(const std::string& arg, const std::string& doc,
                          const bool def);
  /*!
  Get a boolean from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  */
  virtual bool getDocBool(const std::string& arg, const std::string& doc);
  /*!
  Get integer from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param nvals Number of values
  */
  virtual std::vector<int> getDocInts(const std::string& arg,
                                      const std::string& doc, const int nvals);
  /*!
  Get integer from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param defs Default values
  */
  virtual std::vector<int> getDocInts(const std::string& arg,
                                      const std::string& doc,
                                      const std::vector<int>& defs);
  /*!
  Get an floats from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param nval Number of values to look for
  */
  virtual std::vector<float> getDocFloats(const std::string& arg,
                                          const std::string& doc, int nvals);
  /*!
   Add additional parameters

   \param pars List of parameters to add
   */
  virtual void addParams(std::map<std::string, std::string>& pars);

  /*!
   Add additional parameters

   \param pars List of parameters to add = list
   */
  virtual void addParams(std::vector<std::string>& pars);

  /*!
  Get floats from a file

  \param arg Name of the prameter

  \param doc Documentation for parameter

  \param def Default value (if not found in file)
  */
  virtual std::vector<float> getDocFloats(const std::string& arg,
                                          const std::string& doc,
                                          const std::vector<float>& defs);
};  // namespace SEP
}  // namespace SEP

#endif
