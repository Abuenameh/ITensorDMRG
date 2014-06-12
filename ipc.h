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
    cerr << "Write specialization for MPS" << endl;
    stringstream ss(std::ios_base::out);
    mps.write(ss);
    int len = ss.str().length();
    cout << "MPS len: " << len << endl;
    //write(os, len);
    os.write(ss.str().data(), len);
    os.flush();
}

void write(std::ostream& os, IQMPO& mpo) {
    cerr << "Write specialization for MPO" << endl;
    cerr << "Writing MPO" << endl;
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    //write(os, len);
    os.write(ss.str().data(), len);
    os.flush();
    cerr << "Wrote MPO" << endl;
}

template<class T>
void write(std::ostream& os, T& t) {
    cerr << "Write value of type " << type_name<T>() << endl;
    os.write(reinterpret_cast<char*>(&t), sizeof(T));
    os.flush();
}

template<class T>
void write(std::ostream& os, std::vector<T>& v) {
    cerr << "Write vector of type " << type_name<T>() << endl;
    int len = v.size();
    write(os, len);
    for(int i = 0; i < len; i++) {
        write(os, v[i]);
    }
}

/*template<class Tensor>
void write(std::ostream& os, std::vector<MPSt<Tensor> >& v) {
    int len = v.size();
    write(os, len);
    for(int i = 0; i < len; i++) {
        write(os, v[i]);
    }
}*/

void read(std::istream& is, BoseHubbardSiteSet& sites) {
    sites.read(is);
}

void read(std::istream& is, IQMPS& mps) {
    cerr << "Read specialization for MPS" << endl;
    mps.read(is);
}

void read(std::istream& is, IQMPO& mpo) {
    cerr << "Read specialization for MPO" << endl;
    mpo.read(is);
}

template<class T>
void read(std::istream& is, T& t) {
    cerr << "Read value of type " << type_name<T>() << endl;
    is.read(reinterpret_cast<char*>(&t), sizeof(T));
}

template<class T>
void read(std::istream& is, std::vector<T>& v) {
    cerr << "Read vector of type " << type_name<T>() << endl;
    int len;
    read(is, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(is, v[i]);
    }
}

/*template<class Tensor>
void read(std::istream& is, std::vector<MPSt<Tensor> >& v) {
    int len;
    read(is, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(is, v[i]);
    }
}*/

#endif
