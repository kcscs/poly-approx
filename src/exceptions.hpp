#pragma once

// SOURCE: https://stackoverflow.com/questions/24469927/does-c-have-an-equivalent-to-nets-notimplementedexception

#include <stdexcept>

class NotImplemented : public std::logic_error
{
private:

    std::string _text;

    NotImplemented(const char* message, const char* function);

public:

    NotImplemented();

    NotImplemented(const char* message);

    virtual const char *what() const throw();
};
