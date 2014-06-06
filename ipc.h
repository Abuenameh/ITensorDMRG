#ifndef IPC_H
#define IPC_H

#include <sstream>

#include <boost/iostreams/stream.hpp>

#include <ipc_channel.h>
#include <ipc_msg_dispatch.h>
#include <ipc_codec.h>
#include <pipe_unix.h>

#include "BoseHubbardSiteSet.h"

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

void read(std::istream& is, BoseHubbardSiteSet& sites) {
    sites.read(is);
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

typedef ipc::Channel<PipeTransport, ipc::Encoder, ipc::Decoder> PipeChannel;

const int kSitesSend = 1;
const int kRealSend = 2;
const int kRealVectorSend = 3;
const int kIntVectorSend = 4;

class MessageHandler {
//class MessageHandler : public MsgHandlerBaseT {
public:
  bool OnMsgArgCountError(int /*count*/) {
    return false;
  }

  bool OnMsgArgConvertError(int /*code*/) {
    return false;
  }
};

DEFINE_IPC_MSG_CONV(kSitesSend, 1)
{
    IPC_MSG_P1(ipc::ByteArray, ByteArray)
};

class SitesMsgSend : public ipc::MsgOut<PipeChannel>
{
public:
    size_t Send(PipeChannel* ch, BoseHubbardSiteSet& sites) {
        stringstream ss(stringstream::in | stringstream::out | stringstream::binary);
        sites.write(ss);
        return SendMsg(kSitesSend, ch, ipc::ByteArray(ss.str().length(), ss.str().data()));
    }
};

class SitesMsgRecv : public MessageHandler, public ipc::MsgIn<kSitesSend, SitesMsgRecv, PipeChannel>
{
public:
    SitesMsgRecv(BoseHubbardSiteSet& sites) : sites_(sites) {
    }

    virtual size_t OnMsgIn(int msg_id, PipeChannel* ch, const ipc::WireType* const args[], int count) {
        return OnMsgInX(msg_id, ch, args, count);
    }

    size_t OnMsg(PipeChannel*, const ipc::ByteArray& ba) {
        string messagestr;
    messagestr.append(reinterpret_cast<const char*>(ba.buf_), ba.sz_);
    stringstream ss(stringstream::in | stringstream::out | stringstream::binary);
    ss.str(messagestr);
        sites_.read(ss);
        return ipc::OnMsgReady;
    }

private:
    BoseHubbardSiteSet& sites_;
};

DEFINE_IPC_MSG_CONV(kRealSend, 1)
{
    IPC_MSG_P1(ipc::ByteArray, ByteArray)
};

DEFINE_IPC_MSG_CONV(kRealVectorSend, 1)
{
    IPC_MSG_P1(ipc::ByteArray, ByteArray)
};

DEFINE_IPC_MSG_CONV(kIntVectorSend, 1)
{
    IPC_MSG_P1(ipc::ByteArray, ByteArray)
};


#endif
