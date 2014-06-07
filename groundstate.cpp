#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <sstream>
#include <core.h>
#include <hambuilder.h>

#include <boost/multi_array.hpp>
#include <boost/process.hpp>

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
using namespace boost::process;
using namespace boost::process::initializers;

using namespace itensor;

int main(int argc, char **argv)
{
    BoseHubbardSiteSet sites;
    read(cin, sites);
    const int L = sites.N();
    //const int nmax = sites.nmax();
    
    int test;
    //read(cin, test);
    //cerr << "Test = " << test << endl;
    
    IQMPS m(sites);
    //read(cin, m);
    //cerr << m << endl;
    //exit(0);
    
    IQMPO b(sites);
    read(cin, b);
    cerr << b << endl;
    vector<IQMPO> bs(L, IQMPO(sites));
    /*for(int i = 0; i < L; i++) {
        bs.emplace_back(sites);
    }*/
    //read(cin, bs);
    //cerr << bs[1] << endl;

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

            cerr << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << endl;

            BoseHubbardObserver<IQTensor> observer(psi,Opt("EnergyErrgoal",errgoal));

            E0 = dmrg(psi,H,sweeps,observer,Opt("Quiet",quiet));

            psi0 = psi;

            cerr << format("\nGround State Energy = %.10f",E0) << endl;

            Ei = observer.getEnergies();

        } catch(...) {
        }

        time_point<system_clock> end = system_clock::now();
        int runtime = duration_cast<seconds>(end - start).count();

        write(cout, psi0);
        write(cout, E0);
        write(cout, Ei);
        write(cout, runtime);
    }

    return 0;
}
