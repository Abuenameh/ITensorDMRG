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

message_queue *oqueue;

void sigabort(int)
{
    vector<char> buf;
    oqueue->send(buf.data(), 0, 0);
}

int main(int argc, char **argv)
{
    message_queue iq(open_only, argv[1]);
    message_queue oq(open_only, argv[2]);
    
    oqueue = &oq;
    
    signal(SIGINT, sigabort);
    signal(SIGABRT, sigabort);

    std::random_device rd;
    std::mt19937 g(rd());

    BoseHubbardSiteSet sites;
    read(iq, sites);
    const int L = sites.N();

    vector<IQTensor> ns, n2s;
    vector<IQMPO> Cds(L-1, IQMPO(sites));
    for(int i = 1; i <= L; i++) {
        ns.push_back(sites.op("N", i));
        n2s.push_back(sites.op("N2", i));
    }

    read(iq, Cds);

    int nsweeps;
    read(iq, nsweeps);

    Sweeps sweeps(nsweeps);

    vector<int> minm(nsweeps), maxm(nsweeps), niter(nsweeps);
    vector<Real> cutoff(nsweeps), noise(nsweeps);

    read(iq, minm);
    read(iq, maxm);
    read(iq, niter);
    read(iq, cutoff);
    read(iq, noise);
    for(int sw = 1; sw <= nsweeps; sw++) {
        sweeps.setminm(sw, minm[sw-1]);
        sweeps.setmaxm(sw, maxm[sw-1]);
        sweeps.setniter(sw, niter[sw-1]);
        sweeps.setcutoff(sw, cutoff[sw-1]);
        sweeps.setnoise(sw, noise[sw-1]);
    }

    Real errgoal;
    read(iq, errgoal);

    bool quiet;
    read(iq, quiet);
    
    while(true) {

        vector<Real> Js(L), Us(L), mus(L);

        read(iq, Js);
        read(iq, Us);
        read(iq, mus);

        int N;
        read(iq, N);

        time_point<system_clock> start = system_clock::now();

        vector<IQMPS> psis;
        vector<Real> E0s;
        vector<vector<Real> > Eis;

        IQMPS psi0;
        Real E0 = NAN;
        vector<Real> Ei;

            BoseHubbardHamiltonian BH = BoseHubbardHamiltonian(sites);
            BH.J(Js);
            BH.U(Us);
            BH.mu(mus);

            IQMPO H = BH;

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
