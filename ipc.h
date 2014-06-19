#ifndef IPC_H
#define IPC_H

#include <sstream>
#include <cxxabi.h>

#include <zmq.hpp>

#include <boost/iostreams/stream.hpp>

#define BOOST_DATE_TIME_NO_LIB
#include <boost/interprocess/ipc/message_queue.hpp>

#define MAX_MSG_SIZE 100*1024*1024

#include <core.h>
#include "BoseHubbardSiteSet.h"

#undef SP
#undef TCP_NODELAY
#include <nnxx/message>
#include <nnxx/socket>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ostream;
using std::stringstream;

using namespace boost::interprocess;

using namespace zmq;

class connection_monitor : public monitor_t {
public:
    connection_monitor(function<void()> func) : func(func) {}
    
    virtual void on_event_disconnected(const zmq_event_t &event_, const char* addr_) {
        func();
    }
    
private:
    function<void()> func;
};

/*class connection_monitor : public monitor_t {
public:
    connection_monitor(context_t& context, socket_t& socket) : context(context), socket(socket) {}
    
    virtual void on_event_disconnected(const zmq_event_t &event_, const char* addr_) {
        context.close();
    }
    
private:
    context_t& context;
    socket_t& socket;
};*/

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





void write(socket_t& os, BoseHubbardSiteSet& sites) {
    stringstream ss(std::ios_base::out);
    sites.write(ss);
    int len = ss.str().length();
    os.send(ss.str().data(), len);
}

template<class Tensor>
void write(socket_t& os, MPOt<Tensor>& mpo) {
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    os.send(ss.str().data(), len);
}

template<class T>
void write(socket_t& os, T& t) {
    os.send(reinterpret_cast<char*>(&t), sizeof(T));
}

template<class T>
void write(socket_t& os, vector<T>& v) {
    int len = v.size();
    write(os, len);
    for(int i = 0; i < len; i++) {
        write(os, v[i]);
    }
}

void read(socket_t& is, BoseHubbardSiteSet& sites) {
    message_t mesg;
    is.recv(&mesg);
    int len = mesg.size();
    string str = "";
    str.append(reinterpret_cast<char*>(mesg.data()), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
    sites.read(ss);
}

template<class Tensor>
void read(socket_t& is, MPOt<Tensor>& mpo) {
    message_t mesg;
    is.recv(&mesg);
    int len = mesg.size();
    string str = "";
    str.append(reinterpret_cast<char*>(mesg.data()), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
    mpo.read(ss);
}

template<class T>
void read(socket_t& is, T& t) {
    is.recv(reinterpret_cast<char*>(&t), sizeof(T));
}

template<class T>
void read(socket_t& is, vector<T>& v) {
    int len;
    read(is, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(is, v[i]);
    }
}




void write(nnxx::socket& os, BoseHubbardSiteSet& sites) {
    stringstream ss(std::ios_base::out);
    sites.write(ss);
    int len = ss.str().length();
//    cout << "Writing sites len = " << len << endl;
    os.send(ss.str().data(), len, 0);
//    nnxx::message_ostream os;
//    sites.write(os);
//    cout << "Write message length = " << os.msg().size() << endl;
//    s.send(os.msg());
}

void write(nnxx::socket& os, IQMPS& mps) {
    stringstream ss(std::ios_base::out);
    mps.write(ss);
    int len = ss.str().length();
    os.send(ss.str().data(), len, 0);
//    nnxx::message_ostream os;
//    mps.write(os);
//    s.send(os.msg());
}

template<class Tensor>
void write(nnxx::socket& os, MPOt<Tensor>& mpo) {
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    os.send(ss.str().data(), len, 0);
//    nnxx::message_ostream os;
//    mpo.write(os);
//    s.send(os.msg());
}

template<class T>
void write(nnxx::socket& os, T& t) {
//    os.send(static_cast<void*>(&t), sizeof(T), 0);
    os.send(&t, sizeof(T), 0);
}

template<class T>
void write(nnxx::socket& os, vector<T>& v) {
    int len = v.size();
    write(os, len);
    for(int i = 0; i < len; i++) {
        write(os, v[i]);
    }
}

void read(nnxx::socket& is, BoseHubbardSiteSet& sites) {
    nnxx::message msg = is.recv();
    int len = msg.size();
//    cout << "Read sites len = " << len << endl;
    string str = "";
    str.append(static_cast<char*>(msg.data()), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
//    cout << "Sites: " << *((int*)msg.data()) << endl;
    sites.read(ss);
////    nnxx::message_istream is(s.recv());
//    nnxx::message msg = s.recv();
//    cout << "Read message length = " << msg.size() << endl;
//    nnxx::message_istream is(std::move(msg));
//    sites.read(is);
}

template<class Tensor>
void read(nnxx::socket& is, MPOt<Tensor>& mpo) {
    nnxx::message msg = is.recv();
    int len = msg.size();
    string str = "";
    str.append(static_cast<char*>(msg.data()), len);
    stringstream ss(std::ios_base::in);
    ss.str(str);
    mpo.read(ss);
//    nnxx::message_istream is(s.recv());
//    mpo.read(is);
}

template<class T>
void read(nnxx::socket& is, T& t) {
//    is.recv(static_cast<void*>(&t), sizeof(T));
    is.recv(&t, sizeof(T));
}

template<class T>
void read(nnxx::socket& is, vector<T>& v) {
    int len = 0;
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
    mq.send(ss.str().data(), len, 0);
}

template<class Tensor>
void write(message_queue& mq, MPOt<Tensor>& mpo) {
    stringstream ss(std::ios_base::out);
    mpo.write(ss);
    int len = ss.str().length();
    cout << "Writing MPO, len = " << len << endl;
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

void read(message_queue& mq, BoseHubbardSiteSet& sites) {
    unsigned int priority;
    message_queue::size_type len;
    vector<char> buf(MAX_MSG_SIZE);
    mq.receive(buf.data(), MAX_MSG_SIZE, len, priority);
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
    cout << "Read MPO, len = " << len << endl;
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
    mq.receive(&t, sizeof(T), len, priority);
    cout << "Received T " << type_name<T>() << " " << t << endl;
}

template<class T>
void read(message_queue& mq, vector<T>& v) {
    int len;
    cout << "Num msg " << mq.get_num_msg() << endl;
    cout << "About to read vector length" << endl;
    read(mq, len);
    v.resize(len);
    for(int i = 0; i < len; i++) {
        read(mq, v[i]);
    }
}


#endif
