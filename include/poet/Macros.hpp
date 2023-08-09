//  Time-stamp: "Last modified 2023-08-09 12:10:57 delucia"
#ifndef MACROS_H
#define MACROS_H

// Prepend "msg" with name of calling function 
#define MSG(msg)    std::cout << "CPP: " << __func__ << ": " << std::string(msg) << std::endl;

#define MSG_NOENDL(msg)    std::cout << "CPP: " << __func__ << ": " << std::string(msg);

#define ERRMSG(msg) std::cerr << "CPP_ERROR: " << __func__ << ": " << std::string(msg) << std::endl;

#define BOOL_PRINT(bool) std::string(bool ? "ON" : "OFF")

#endif // MACROS_H
