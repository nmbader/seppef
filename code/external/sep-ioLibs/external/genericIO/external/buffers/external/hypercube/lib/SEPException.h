#ifndef sepException_H
#define sepException_H
#include <string>
/*!
A simple exception class.
Tied in with pybind to allow c++ exceptions to be thrown in python*/
class SEPException : public std::exception {
 public:
  explicit SEPException(std::string m)
      : message{m} {}  ///< Create exception from a string
  virtual const char* what() const noexcept override {
    return message.c_str();
  }  ///< Return the exception char
  std::string getMessage() const {
    return message;
  }  ///< Return exception as string
 private:
  std::string message = "";  ///< Storage for exception text
};
#endif
