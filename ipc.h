#ifndef IPC_H
#define IPC_H

#include <sstream>
#include <cxxabi.h>

#include <boost/iostreams/stream.hpp>

#include "BoseHubbardSiteSet.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::stringstream;

template<typename T>
std::string type_name()
{
    int status;
    std::string tname = typeid(T).name();
    char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if(status == 0) {
        tname = demangled_name;
        std::free(demangled_name);
    }   
    return tname;
}

void write(std::ostream& os, BoseHubbardSiteSet& sites) {
    sites.write(os);
    os.flush();
}

void write(std::ostream& os, IQMPS& mps) {
    stringstream ss(std::ios_base::out);
    mps.write(ss);
    int len = ss.str().length();
    os.write(ss.str().data(), len);
    os.flush();
}

void write(std::ostream& os, IQMPO& mpo) {
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    os.write(ss.str().data(), len);
    os.flush();
}

template<class T>
void write(std::ostream& os, T& t) {
    os.write(reinterpret_cast<char*>(&t), sizeof(T));
    os.flush();
}

template<class T>
void write(std::ostream& os, std::vector<T>& v) {
    int len = v.size();
    write(os, len);
    for(int i = 0; i < len; i++) {
        write(os, v[i]);
    }
}

void read(std::istream& is, BoseHubbardSiteSet& sites) {
    sites.read(is);
}

void read(std::istream& is, IQMPS& mps) {
    mps.read(is);
}

void read(std::istream& is, IQMPO& mpo) {
    mpo.read(is);
}

template<class T>
void read(std::istream& is, T& t) {
    is.read(reinterpret_cast<char*>(&t), sizeof(T));
}

template<class T>
void read(std::istream& is, std::vector<T>& v) {
    int len;
    read(is, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(is, v[i]);
    }
}

#endif
