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

using std::ref;
using std::stoi;
using std::stod;
using std::cin;
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
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

using namespace itensor;

stream<file_descriptor_sink> *abortos;

void sigabort(int)
{
    bool abort = true;
    write(*abortos, abort);
    exit(1);
}

int main(int argc, char **argv)
{
    stream<file_descriptor_sink> os(42, never_close_handle);
    abortos = &os;
    signal(SIGABRT, sigabort);
    signal(SIGINT, sigabort);

    BoseHubbardSiteSet sites;
    read(cin, sites);
    const int L = sites.N();

    vector<IQTensor> ns, n2s;
    vector<vector<IQMPO> > bdbs(L, vector<IQMPO>(L));
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
    read(cin, Cds);

    int nsweeps;
    read(cin, nsweeps);

    Sweeps sweeps(nsweeps);

    vector<int> minm(nsweeps), maxm(nsweeps), niter(nsweeps);
    vector<Real> cutoff(nsweeps), noise(nsweeps);

    read(cin, minm);
    read(cin, maxm);
    read(cin, niter);
    read(cin, cutoff);
    read(cin, noise);
    for(int sw = 1; sw <= nsweeps; sw++) {
        sweeps.setminm(sw, minm[sw-1]);
        sweeps.setmaxm(sw, maxm[sw-1]);
        sweeps.setniter(sw, niter[sw-1]);
        sweeps.setcutoff(sw, cutoff[sw-1]);
        sweeps.setnoise(sw, noise[sw-1]);
    }

    Real errgoal;
    read(cin, errgoal);

    bool quiet;
    read(cin, quiet);

    while(true) {

        vector<Real> Js(L), Us(L), mus(L);

        read(cin, Js);
        read(cin, Us);
        read(cin, mus);

        int N;
        read(cin, N);

        time_point<system_clock> start = system_clock::now();

        vector<IQMPS> psis;
        vector<Real> E0s;
        vector<vector<Real> > Eis;

        IQMPS psi0;
        Real E0 = NAN;
        vector<Real> Ei;

        try {

            BoseHubbardHamiltonian BH = BoseHubbardHamiltonian(sites);
            BH.J(Js);
            BH.U(Us);
            BH.mu(mus);

            IQMPO H = BH;

            for(int eig = 0; eig < 1; eig++) {

                InitState initState(sites);
                int n0 = N / L;
                for(int i = 1; i <= L; ++i) {
                    if(i <= N % L)
                        initState.set(L+1-i,std::to_string(n0+1));
                    else
                        initState.set(L+1-i,std::to_string(n0));
                }
                IQMPS psi(initState);

                cerr << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << endl;

                BoseHubbardObserver<IQTensor> observer(psi,Opt("EnergyErrgoal",errgoal));

                E0 = dmrg(psi,H,psis,sweeps,observer,Opt("Quiet",quiet));

                psi0 = psi;

                cerr << format("\nGround State Energy = %.10f",E0) << endl;

                Ei = observer.getEnergies();

                psis.push_back(psi0);
                E0s.push_back(E0);
                Eis.push_back(Ei);

            }

        } catch(...) {
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

        write(cout, E0);
        write(cout, Ei);
        write(cout, n);
        write(cout, n2);
        write(cout, C);
        write(cout, runtime);

        bool abort = false;
        write(*abortos, abort);
    }

    return 0;
}
