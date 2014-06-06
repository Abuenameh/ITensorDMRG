#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <sstream>
#include <core.h>
#include <hambuilder.h>

#include <boost/multi_array.hpp>
#include <boost/process.hpp>

#include <zmq.hpp>

#include "ipc.h"

#include "concurrent_queue.h"
#include "ThreadPool.h"
#include "mathematica.h"

#include "BoseHubbardSiteSet.h"
#include "BoseHubbardHamiltonian.h"
#include "BoseHubbardObserver.h"

using std::ref;
using std::stoi;
using std::stod;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::queue;
using std::to_string;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::atomic_bool;
using std::chrono::time_point;
using std::chrono::duration;
using std::chrono::system_clock;
using std::chrono::seconds;
using std::chrono::duration_cast;

using boost::extents;
using boost::multi_array;
using namespace boost::process;
using namespace boost::process::initializers;

using namespace itensor;

using namespace zmq;

#define ZMQ_MAX_SIZE 10*1024*1024

stringstream get_message_stream(socket_t& socket)
{
    message_t message(ZMQ_MAX_SIZE);
    string messagestr;

    int len = socket.recv(&message);
    messagestr.append(reinterpret_cast<char*>(message.data()), len);

    stringstream ss(stringstream::in | stringstream::out | stringstream::binary);
    ss.str(messagestr);
    return ss;
}

template<class T>
void get_message(socket_t& socket, T& v)
{
    stringstream ss = get_message_stream(socket);
    ss.read(reinterpret_cast<char*>(&v), sizeof(T));
}

template<class T>
void send_message(socket_t& socket, T& v)
{
    message_t message(ZMQ_MAX_SIZE);
    socket.send(&v, sizeof(T));
}

int main(int argc, char **argv)
{
    BoseHubbardSiteSet qwe;
    cerr << "About to read" << endl;
    //qwe.read(std::cin);
    read(std::cin, qwe);
    cerr << qwe << endl;
    cerr << "Read" << endl;
    //int d = std::cin.get();
    //cerr << d << endl;
    int a, b, c;
    //std::cin >> a >> a;
    read(std::cin, a);
    vector<double> v;
    read(std::cin, v);
    for(int i = 0; i < v.size(); i++) {
        cerr << v[i] << " ";
    }
    cerr << endl;
    //cerr << a << endl;
    //std::cin >> a >> b >> c;
    //cerr << a << endl << b << endl << c << endl;
    return 0;
    
    int port = stoi(argv[1]);

    context_t context;
    socket_t socket(context, ZMQ_PULL), socketout(context, ZMQ_PUSH);

    string addr = "tcp://localhost:" + to_string(port);
    socket.bind(addr.c_str());
    string addrout = "tcp://localhost:" + to_string(port+1);
    socketout.connect(addrout.c_str());

    message_t message(ZMQ_MAX_SIZE);
    string messagestr;

    stringstream ss;

        ss = get_message_stream(socket);
        BoseHubbardSiteSet sites;
        sites.read(ss);
        cout << sites << endl;

        const int L = sites.N();
        //const int nmax = sites.nmax();

    while(true) {
        /*stringstream ss(stringstream::in | stringstream::out | stringstream::binary);

        int len = socket.recv(&message);
        messagestr.clear();
        messagestr.append(reinterpret_cast<char*>(message.data()), len);
        ss.str(messagestr);*/

        /*ss = get_message_stream(socket);
        BoseHubbardSiteSet sites;
        sites.read(ss);

        const int L = sites.N();
        //const int nmax = sites.nmax();*/

        vector<Real> Js(L), Us(L), mus(L);

        //ss = get_message_stream(socket);
        for(int i = 0; i < L; i++) {
            //Real J;
            //ss.read(reinterpret_cast<char*>(&Js[i]), sizeof(Real));
            get_message(socket, Js[i]);
        }
        //ss = get_message_stream(socket);
        for(int i = 0; i < L; i++) {
            //Real U;
            //ss.read(reinterpret_cast<char*>(&Us[i]), sizeof(Real));
            get_message(socket, Us[i]);
        }
        //ss = get_message_stream(socket);
        for(int i = 0; i < L; i++) {
            //Real mu;
            //ss.read(reinterpret_cast<char*>(&mus[i]), sizeof(Real));
            get_message(socket, mus[i]);
        }
        
        int nsweeps;
        get_message(socket, nsweeps);
        
        Sweeps sweeps(nsweeps);

        //vector<int> minm(nsweeps), maxm(nsweeps), niter(nsweeps);
        //vector<Real> cutoff(nsweeps), noise(nsweeps);
        int minm, maxm, niter;
        Real cutoff, noise;
        for(int sw = 0; sw < nsweeps; sw++) {
            //get_message(socket, minm[sw]);
            get_message(socket, minm);
            sweeps.setminm(sw, minm);
        }
        for(int sw = 0; sw < nsweeps; sw++) {
            //get_message(socket, maxm[sw]);
            get_message(socket, maxm);
            sweeps.setmaxm(sw, maxm);
        }
        for(int sw = 0; sw < nsweeps; sw++) {
            //get_message(socket, niter[sw]);
            get_message(socket, niter);
            sweeps.setniter(sw, niter);
        }
        for(int sw = 0; sw < nsweeps; sw++) {
            //get_message(socket, cutoff[sw]);
            get_message(socket, cutoff);
            sweeps.setcutoff(sw, cutoff);
        }
        for(int sw = 0; sw < nsweeps; sw++) {
            //get_message(socket, noise[sw]);
            get_message(socket, noise);
            sweeps.setnoise(sw, noise);
        }
        
        int N;
        get_message(socket, N);
        //ss = get_message(socket);
        //ss.read(reinterpret_cast<char*>(&N), sizeof(int));
        
        bool quiet;
        get_message(socket, quiet);
        
        double errgoal;
        get_message(socket, errgoal);

        time_point<system_clock> start = system_clock::now();

        IQMPS psi0;
        Real E0 = NAN;
        vector<Real> Ei;

        try {

            BoseHubbardHamiltonian BH = BoseHubbardHamiltonian(sites);
            BH.J(Js);
            BH.U(Us);
            BH.mu(mus);

            IQMPO H = BH;

            InitState initState(sites);
            int n0 = N / L;
            for(int i = 1; i <= L; ++i) {
                if(i <= N % L)
                    initState.set(i,std::to_string(n0+1));
                else
                    initState.set(i,std::to_string(n0));
            }
            IQMPS psi(initState);

            cout << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << endl;

            BoseHubbardObserver<IQTensor> observer(psi,Opt("EnergyErrgoal",errgoal));

            E0 = dmrg(psi,H,sweeps,observer,Opt("Quiet",quiet));

            psi0 = psi;

            cout << format("\nGround State Energy = %.10f",E0) << endl;

            Ei = observer.getEnergies();

        } catch(...) {
        }

        time_point<system_clock> end = system_clock::now();
        int runtime = duration_cast<seconds>(end - start).count();

        stringstream ssout(stringstream::in | stringstream::out | stringstream::binary);
        psi0.write(ssout);
        socketout.send(ssout.str().data(), ssout.str().size());

        send_message(socketout, E0);
        
        int Eilen = Ei.size();
        send_message(socketout, Eilen);
        socketout.send(Ei.data(), Eilen*sizeof(Real));
        
        send_message(socketout, runtime);
    }

    return 0;
}
