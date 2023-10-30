#ifndef gcp_buffers_reg_file_h
#define gcp_buffers_reg_file_h 1
#include <stdbool.h>
#include <string>
#include "basicIO.h"
#include "buffers.h"
#include "buffersRegFile.h"
#include "google/cloud/status_or.h"
#include "google/cloud/storage/client.h"
#include "json.h"
namespace SEP {
/*!
Object for handling regular files broken into blocks and written to objects on
GCP
*/
class gcpBuffersRegFile : public buffersRegFile {
 public:
  // sepRegFile::sepRegFile(const std::string tag,usage_code usage){
  /*!
      Initialize the default GCP boject
  */
  gcpBuffersRegFile() { ; }
  /*!
    Initialize GCP regular file
      \param arg  JSON arguments
      \param usage Usage for file
      \param tag Tag associated with GCP file
      \param progName Program name
      \param ndimMax Minimum number of dimensions
*/
  gcpBuffersRegFile(const Json::Value &arg, const SEP::usage_code usage,
                    const std::string &tag, const std::string &progName,
                    const int ndimMax);
  /*!
    Setup GCP environment

    \param arg  Arguments (JSON) which describes file
    \param tag Tag associated with file
    */
  void setupGCP(const Json::Value &arg, const std::string &tag);

  /*!
   Remove Description and directory
   */
  virtual void removeDescDir() override;

  /*!
    Close file
    */

  virtual void close();
  /*!
      Create buffers
  */
  void createBuffers();

  std::string getEnvVar(std::string const &key,
                        std::string const &defaultV) const {
    char *val = getenv(key.c_str());
    return val == NULL ? defaultV : std::string(val);
  }

 protected:
  google::cloud::v0::StatusOr<google::cloud::storage::Client> _client;

 private:
  /*!
Setup GCP
*/
  void setupGCP();
  std::string _bucket;  ///< Name of the bucker
  std::string _dir;     ///< Directory

};  // namespace SEP

}  // namespace SEP
#endif
