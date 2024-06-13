#ifndef MACROS_H
#define MACROS_H

#include <iostream>
#include <string>

// Prepend "msg" with name of calling function 
#define MSG(msg)    std::cout << "CPP: " << __func__ << ": " << std::string(msg) << std::endl;

#define MSG_NOENDL(msg)    std::cout << "CPP: " << __func__ << ": " << std::string(msg);

#define ERRMSG(msg) std::cerr << "CPP_ERROR: " << __func__ << ": " << std::string(msg) << std::endl;

#define BOOL_PRINT(bool) std::string(bool ? "ON" : "OFF")

#endif // MACROS_H
