#include "exceptions.hpp"

NotImplemented::NotImplemented(const char *message, const char *function)
    : std::logic_error("Not Implemented") {
  _text = message;
  _text += " : ";
  _text += function;
};

NotImplemented::NotImplemented()
    : NotImplemented("Not Implememented", __FUNCTION__) {}

NotImplemented::NotImplemented(const char *message)
    : NotImplemented(message, __FUNCTION__) {}

const char *NotImplemented::what() const throw() { return _text.c_str(); }
