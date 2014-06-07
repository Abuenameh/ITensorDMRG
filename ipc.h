#ifndef IPC_H
#define IPC_H

#include <sstream>

#include <boost/iostreams/stream.hpp>

#include "BoseHubbardSiteSet.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::stringstream;

void write(std::ostream& os, BoseHubbardSiteSet& sites) {
    sites.write(os);
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

void write(std::ostream& os, IQMPS& mps) {
    stringstream ss(std::ios_base::out);
    mps.write(ss);
    int len = ss.str().length();
    write(os, len);
    os.write(ss.str().data(), len);
    os.flush();
}

void write(std::ostream& os, IQMPO& mpo) {
    cerr << "Writing MPO" << endl;
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    write(os, len);
    os.write(ss.str().data(), len);
    os.flush();
    cerr << "Wrote MPO" << endl;
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
