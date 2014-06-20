#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <sstream>
#include <core.h>
#include <hambuilder.h>

#include <boost/multi_array.hpp>
#include <boost/process.hpp>
#include <boost/iostreams/stream.hpp>

#include "ipc.h"

#include "concurrent_queue.h"
#include "ThreadPool.h"
#include "mathematica.h"

#include "BoseHubbardSiteSet.h"
#include "BoseHubbardHamiltonian.h"
#include "BoseHubbardObserver.h"

#include <zmq.hpp>

using std::ref;
using std::stoi;
using std::stod;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::string;
using std::vector;
using std::queue;
using std::thread;
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
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

using namespace itensor;

using namespace zmq;

#ifndef FST

stream<file_descriptor_sink> *abortos;

void sigabort(int)
{
    bool abort = true;
    write(*abortos, abort);
    exit(1);
}

#else

//char *oqueue_name;
//char *iqueue_name;
message_queue *oqueue;

void sigabort(int)
{
    cout << endl << "------------------Aborting------------------" << endl;
    vector<char> buf;
    oqueue->send(buf.data(), 0, 0);
}
#endif

int main(int argc, char **argv)
{
    /*string base = "tcp://127.0.0.1:";
    string inport = base + argv[3];
    string outport = base + argv[4];*/
    
    /*context_t context(1);
    socket_t is(context, ZMQ_PULL);
    socket_t os(context, ZMQ_PUSH);
    
    is.bind(inport.c_str());
    os.connect(outport.c_str());*/
    
    /*nnxx::socket is(nnxx::SP, nnxx::PULL);
    nnxx::socket os(nnxx::SP, nnxx::PUSH);
    
    is.bind(inport);
    os.connect(outport);
    
    cout << "Created sockets" << endl;*/
    
//    string queue_name = argv[1];
    message_queue iq(open_only, argv[1]);
    message_queue oq(open_only, argv[2]);
    
//    iqueue_name = argv[1];
//    oqueue_name = argv[2];
    oqueue = &oq;
    
    cout << "Opened message queues" << endl;
    
    
#ifndef FST
    stream<file_descriptor_sink> mq(42, never_close_handle);
    abortos = &mq;
    signal(SIGABRT, sigabort);
    signal(SIGINT, sigabort);
#else
    signal(SIGABRT, sigabort);
#endif

    std::random_device rd;
    std::mt19937 g(rd());

    BoseHubbardSiteSet sites;
    read(iq, sites);
    const int L = sites.N();
    cout << "Read sites" << endl << flush;

    vector<IQTensor> ns, n2s;
    vector<IQMPO> Cds(L-1, IQMPO(sites));
    for(int i = 1; i <= L; i++) {
        ns.push_back(sites.op("N", i));
        n2s.push_back(sites.op("N2", i));
    }

    /*for(int d = 1; d < L; d++) {
        IQMPO Cd = HamBuilder<IQTensor>(sites, "Bdag", 1, "B", 1+d);
        for(int i = 2; i <= L-d; i++) {
            Cd.plusEq(HamBuilder<IQTensor>(sites, "Bdag", i, "B", i+d));
        }
        Cd *= 1./(L-d);
        Cds.push_back(Cd);
    }*/
    cout << "About to read Cds" << endl << flush;
    read(iq, Cds);
    cout << "Read Cds" << endl << flush;

    int nsweeps;
    read(iq, nsweeps);
    cout << "Read nsweeps" << endl << flush;

    Sweeps sweeps(nsweeps);

    vector<int> minm(nsweeps), maxm(nsweeps), niter(nsweeps);
    vector<Real> cutoff(nsweeps), noise(nsweeps);

    read(iq, minm);
    read(iq, maxm);
    read(iq, niter);
    read(iq, cutoff);
    read(iq, noise);
    cout << "Read sweeps" << endl << flush;
    for(int sw = 1; sw <= nsweeps; sw++) {
        sweeps.setminm(sw, minm[sw-1]);
        sweeps.setmaxm(sw, maxm[sw-1]);
        sweeps.setniter(sw, niter[sw-1]);
        sweeps.setcutoff(sw, cutoff[sw-1]);
        sweeps.setnoise(sw, noise[sw-1]);
    }

    Real errgoal;
    read(iq, errgoal);
    cout << "Read errgoal" << endl << flush;

    bool quiet;
    read(iq, quiet);
    
    cout << "Read setup" << endl << flush;

//    int wait = 1;
//    write(oq, wait);
//    abort();
    
    while(true) {

//        cout << "Creating parameter vectors " + inport << endl << flush;
        vector<Real> Js(L), Us(L), mus(L);
//        cout << "Created parameter vectors " + inport << endl << flush;

        try{
        read(iq, Js);
        } catch(...) {cout << "----------------------Exception----------------------" << endl << flush;}
//        cout << "Read J " + inport << endl;
        read(iq, Us);
//        cout << "Read U " + inport << endl;
        read(iq, mus);
//        cout << "Read parameters " + inport << endl;

        int N;
        read(iq, N);
        cout << "Read particle number" << endl;

        time_point<system_clock> start = system_clock::now();

        vector<IQMPS> psis;
        vector<Real> E0s;
        vector<vector<Real> > Eis;

        IQMPS psi0;
        Real E0 = NAN;
        vector<Real> Ei;

//        try {

            BoseHubbardHamiltonian BH = BoseHubbardHamiltonian(sites);
            BH.J(Js);
            BH.U(Us);
            BH.mu(mus);
            cout << "Created Hamiltonian" << endl;

            IQMPO H = BH;
            cout << "Converted Hamiltonian" << endl;

            for(int eig = 0; eig < 1; eig++) {

                std::vector<int> ind;
                for(int i = 1; i <= L; i++)
                    ind.push_back(i);
                std::shuffle(ind.begin(), ind.end(), g);
                
                InitState initState(sites);
                int n0 = N / L;
                for(int i = 1; i <= L; ++i) {
                    if(i <= N % L)
                        initState.set(ind[i-1],std::to_string(n0+1));
                    else
                        initState.set(ind[i-1],std::to_string(n0));
                }
                IQMPS psi(initState);

                cerr << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << endl;

                BoseHubbardObserver<IQTensor> observer(psi,Opt("EnergyErrgoal",errgoal));

                E0 = dmrg(psi,H,/*psis,*/sweeps,observer,Opt("Quiet",quiet));

                psi0 = psi;

                cerr << format("\nGround State Energy = %.10f",E0) << endl;

                Ei = observer.getEnergies();

                psis.push_back(psi0);
                E0s.push_back(E0);
                Eis.push_back(Ei);

            }

//        } catch(...) {
//        }

        auto minE0 = min_element(E0s.begin(), E0s.end());
        int pos = distance(E0s.begin(), minE0);
        E0 = E0s[pos];
        psi0 = psis[pos];
        Ei = Eis[pos];

        vector<Real> C;
        for(int i = 0; i < L-1; i++) {
            C.push_back(psiHphi(psi0, Cds[i], psi0));
        }

        vector<Real> n(L), n2(L);
        for(int i = 0; i < L; i++) {
            psi0.position(i+1);
            n[i] = Dot(conj(primed(psi0.A(i+1),Site)),ns[i]*psi0.A(i+1));
            n2[i] = Dot(conj(primed(psi0.A(i+1),Site)),n2s[i]*psi0.A(i+1));
        }

        time_point<system_clock> end = system_clock::now();
        int runtime = duration_cast<seconds>(end - start).count();

        write(oq, E0);
        write(oq, Ei);
        write(oq, n);
        write(oq, n2);
        write(oq, C);
        write(oq, runtime);
    }

    return 0;
}
