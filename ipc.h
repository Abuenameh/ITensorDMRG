#ifndef IPC_H
#define IPC_H

#include <stdexcept>
#include <sstream>
#include <cxxabi.h>

#include <zmq.hpp>

#include <boost/iostreams/stream.hpp>

#define BOOST_DATE_TIME_NO_LIB
#include <boost/interprocess/ipc/message_queue.hpp>

#define NUM_MSG 10
//#define MAX_MSG_SIZE 10*1024*1024
#define MAX_MSG_SIZE 100*1024

//#define SETUP_NUM_MSG L+5*nsweeps+4
#define SETUP_NUM_MSG L
#define SETUP_MAX_MSG_SIZE 10*1024*1024

#define CDS_MAX_MSG_SIZE 10*1024*1024

#include <core.h>
#include "BoseHubbardSiteSet.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ostream;
using std::stringstream;
using std::runtime_error;

using namespace boost::interprocess;

using namespace zmq;

template<typename T>
string type_name()
{
    int status;
    string tname = typeid(T).name();
    char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if(status == 0) {
        tname = demangled_name;
        std::free(demangled_name);
    }   
    return tname;
}

void write(ostream& os, BoseHubbardSiteSet& sites) {
    sites.write(os);
    os.flush();
}

void write(ostream& os, IQMPS& mps) {
    stringstream ss(std::ios_base::out);
    mps.write(ss);
    int len = ss.str().length();
    os.write(ss.str().data(), len);
    os.flush();
}

void write(ostream& os, IQMPO& mpo) {
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    os.write(ss.str().data(), len);
    os.flush();
}

template<class T>
void write(ostream& os, T& t) {
    os.write(reinterpret_cast<char*>(&t), sizeof(T));
    os.flush();
}

template<class T>
void write(ostream& os, vector<T>& v) {
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
void read(std::istream& is, vector<T>& v) {
    int len;
    read(is, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(is, v[i]);
    }
}

void write(message_queue& mq, BoseHubbardSiteSet& sites) {
    stringstream ss(std::ios_base::out);
    sites.write(ss);
    int len = ss.str().length();
    cout << "BH len = " << len << endl;
    mq.send(ss.str().data(), len, 0);
}

template<class Tensor>
void write(message_queue& mq, MPOt<Tensor>& mpo) {
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    cout << "Write MPO len = " << len << endl;
    mq.send(ss.str().data(), len, 0);
}

template<class T>
void write(message_queue& mq, T& t) {
    mq.send(&t, sizeof(T), 0);
}

template<class T>
void write(message_queue& mq, vector<T>& v) {
    int len = v.size();
    write(mq, len);
    for(int i = 0; i < len; i++) {
        write(mq, v[i]);
    }
}

inline void check_termination(int len) {
    if (len == 0) {
        throw runtime_error("Connection terminated");
    }
}

void read(message_queue& mq, BoseHubbardSiteSet& sites) {
    unsigned int priority;
    message_queue::size_type len;
    vector<char> buf(MAX_MSG_SIZE);
    mq.receive(buf.data(), MAX_MSG_SIZE, len, priority);
    check_termination(len);
    string str = "";
    str.append(buf.data(), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
    sites.read(ss);
}

template<class Tensor>
void read(message_queue& mq, MPOt<Tensor>& mpo) {
    unsigned int priority;
    message_queue::size_type len;
    vector<char> buf(MAX_MSG_SIZE);
    mq.receive(buf.data(), MAX_MSG_SIZE, len, priority);
    check_termination(len);
    string str = "";
    str.append(buf.data(), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
    mpo.read(ss);
}

template<class Tensor>
void readCd(message_queue& mq, MPOt<Tensor>& mpo) {
    unsigned int priority;
    message_queue::size_type len;
    vector<char> buf(SETUP_MAX_MSG_SIZE);
    mq.receive(buf.data(), SETUP_MAX_MSG_SIZE, len, priority);
    cout << "Read MPO len = " << len << endl;
    check_termination(len);
    string str = "";
    str.append(buf.data(), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
    mpo.read(ss);
}

template<class T>
void read(message_queue& mq, T& t) {
    unsigned int priority;
    message_queue::size_type len;
    vector<char> buf(MAX_MSG_SIZE);
    mq.receive(buf.data(), MAX_MSG_SIZE, len, priority);
    check_termination(len);
    t = *reinterpret_cast<T*>(buf.data());
}

template<class T>
void readCd(message_queue& mq, T& t) {
    unsigned int priority;
    message_queue::size_type len;
    vector<char> buf(SETUP_MAX_MSG_SIZE);
    mq.receive(buf.data(), SETUP_MAX_MSG_SIZE, len, priority);
    check_termination(len);
    t = *reinterpret_cast<T*>(buf.data());
}

template<class T>
void read(message_queue& mq, vector<T>& v) {
    int len;
    read(mq, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(mq, v[i]);
    }
}

template<class Tensor>
void readCds(message_queue& mq, vector<MPOt<Tensor> >& v) {
    int len;
    readCd(mq, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        readCd(mq, v[i]);
    }
}


#endif
