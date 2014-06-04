#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <core.h>
#include <hambuilder.h>

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
using std::ofstream;
using std::ifstream;
using std::chrono::time_point;
using std::chrono::duration;
using std::chrono::system_clock;

using namespace itensor;

struct Results {
    int it;
    int iN;
    vector<Real> ts;
    vector<Real> Us;
    vector<Real> mus;
    Real E0;
    vector<Real> Ei;
    Vector n;
    Vector n2;
    Matrix C;
    int runtime;
};

void groundstate(concurrent_queue<Results>& resq, Sweeps& sweeps, Real errgoal, bool quiet, BoseHubbardSiteSet& sites, vector<vector<IQMPO> >& COp, int it, int iN, vector<Real> ts, vector<Real> Us, vector<Real> mus, int N)
{
    time_point<system_clock> start = system_clock::now();

    int L = sites.N();

    BoseHubbardHamiltonian BH = BoseHubbardHamiltonian(sites);
    BH.t(ts);
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

    std::cout << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << std::endl;

    BoseHubbardObserver<IQTensor> observer(psi,Opt("EnergyErrgoal",errgoal));

    Real E0 = dmrg(psi,H,sweeps,observer,Opt("Quiet",quiet));

    std::cout << format("\nGround State Energy = %.10f",E0) << std::endl;

    std::vector<Real> Ei = observer.getEnergies();

    Vector n(L), n2(L);
    Matrix C(L, L);
    for(int j = 1; j <= L; ++j) {
        psi.position(j);
        n(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("N",j)*psi.A(j));
        n2(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("N2",j)*psi.A(j));
        for(int k = 1; k <= L; ++k) {
            C(j, k) = psiHphi(psi,COp[j-1][k-1],psi);
        }
    }

    time_point<system_clock> end = system_clock::now();
    duration<double> runtime = end - start;

    Results res;
    res.it = it;
    res.iN = iN;
    res.ts = ts;
    res.Us = Us;
    res.mus = mus;
    res.E0 = E0;
    res.Ei = Ei;
    res.n = n;
    res.n2 = n2;
    res.C = C;
}

vector<Real> linspace(Real min, Real max, int n) {
    vector<Real> result(n);
    double step = (max-min) / (n - 1);
    for(int i = 0; i < n-1; ++i) {
        result[i] = min + i*step;
    }
    result[n-1] = max;
}

int main2(int argc, char **argv)
{
    if(argc < 10) {
        cerr << "Insufficient number of arguments" << endl;
        return 1;
    }
    
    time_point<system_clock> start = system_clock::now();

    int resi = stoi(argv[1]);
    double ximax = stod(argv[2]);
    int seed = stoi(argv[3]);
    
    Real tmin = stod(argv[4]);
    Real tmax = stod(argv[5]);
    int nt = stoi(argv[6]);
    vector<Real> tv = linspace(tmin, tmax, nt);
    
    int Nmin = stoi(argv[7]);
    int Nmax = stoi(argv[8]);
    int nN = Nmax - Nmin + 1;
    vector<int> Nv(nN);
    for(int i = 0; i < nN; ++i) {
        Nv[i] = Nmin + i;
    }
    
    int numthreads = stoi(argv[9]);

    const int L = 10;
    const int nmax = 7;

    //const int nsweeps = 5;
    const Real errgoal = -1;
    const bool quiet = true;

    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;


    BoseHubbardSiteSet sites(L, nmax);

    std::vector<IQIndex> links(L+1);
    for(int l = 0; l <= L; ++l) {
        std::vector<IndexQN> indices;
        for(int n = 0; n <= nmax; n++) {
            indices.push_back(IndexQN(Index(nameint("n = ", n) + nameint(" for site ", l),1), QN(0,n,n%2)));
        }
        links.at(l) = IQIndex(nameint("BoseHubbard site=",l),indices);
    }

    vector<vector<IQMPO> > COp;

    for(int i = 1; i <= L; ++i) {
        vector<IQMPO> Ci;
        for(int j = 1; j <= L; ++j) {
            IQMPO Cij(sites);
            for(int n = 1; n <= L; ++n) {
                IQTensor& W = Cij.Anc(n);
                IQIndex row = dag(links[n-1]), col = links[n];

                W = IQTensor(dag(sites.si(n)),sites.siP(n),row,col);

                if(i != j && (n == i || n == j)) {
                    W += sites.op("Bdag",n) * row(1) * col(2);
                    W += sites.op("B",n) * row(2) * col(4);
                } else {
                    for(int k = 1; k <= 4; ++k) {
                        W += sites.op("Id",n) * row(k) * col(k);
                    }
                }
            }
            Cij.Anc(1) *= IQTensor(links.at(0)(1));
            Cij.Anc(L) *= IQTensor(dag(links.at(L)(4)));
            Ci.push_back(Cij);
        }
        COp.push_back(Ci);
    }
    
    ThreadPool pool(numthreads);
    concurrent_queue<Results> resq;
    
    vector<Real> Us(L, 1);
    vector<Real> mus(L, 0);
    
    for(int it = 0; it < nt; ++it) {
        for(int iN = 0; iN < nN; ++iN) {
            vector<Real> ts(L, tv[it]);
            //pool.enqueue([=,&resq,&sweeps,&sites,&COp] { groundstate(resq, sweeps, errgoal, quiet, sites, COp, it, iN, ts, Us, mus, Nv[iN]); });
            pool.enqueue(bind(groundstate, ref(resq), ref(sweeps), errgoal, quiet, ref(sites), ref(COp), it, iN, ts, Us, mus, Nv[iN]));
        }
    }
    
    Results res;
    while(true) {
        resq.wait_and_pop(res);
        cout << math(res.n) << endl;
    }

    time_point<system_clock> end = system_clock::now();
    std::chrono::seconds runtime = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "runtime = " << runtime.count() << endl;

    return 0;
}

