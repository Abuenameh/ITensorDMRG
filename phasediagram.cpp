#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <core.h>
#include <hambuilder.h>

#include <boost/multi_array.hpp>
#include <boost/process.hpp>

#include <zmq.hpp>

#include "ipc.h"
#include "concurrent_queue.h"
#include "ThreadPool.h"
#include "mathematica.h"
#include "ProcessPool.h"

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

struct Results {
    int ix;
    int iN;
    vector<Real> xs;
    vector<Real> Us;
    vector<Real> mus;
    Real E0;
    vector<Real> Ei;
    vector<Real> n;
    vector<Real> n2;
    vector<vector<Real> > C;
    int runtime;
};

string seconds_to_string(int s)
{
    int m = s / 60;
    s = s % 60;
    int h = m / 60;
    m = m % 60;
    int d = h / 24;
    h = h % 24;

    if(d > 0) {
        return format("%d d %d:%02d:%02d", d, h, m, s);
    } else {
        return format("%d:%02d:%02d", h, m, s);
    }
}


vector<Real> linspace(Real min, Real max, int n)
{
    vector<Real> result(n);
    double step = (max-min) / (n - 1);
    for(int i = 0; i < n-1; ++i) {
        result[i] = min + i*step;
    }
    result[n-1] = max;
    return result;
}

atomic_bool interrupted;

void sigint(int)
{
    interrupted = true;
}

int main(int argc, char **argv)
{
    if(argc < 10) {
        cerr << "Insufficient number of arguments" << endl;
        return 1;
    }

    time_point<system_clock> start = system_clock::now();

    int resi = stoi(argv[1]);

    double ximax = stod(argv[2]);
    int seed = stoi(argv[3]);

    const Real xmin = stod(argv[4]);
    const Real xmax = stod(argv[5]);
    const int nx = stoi(argv[6]);
    vector<Real> xv = linspace(xmin, xmax, nx);

    const int Nmin = stoi(argv[7]);
    const int Nmax = stoi(argv[8]);
    const int nN = Nmax - Nmin + 1;
    vector<int> Nv(nN);
    for(int i = 0; i < nN; ++i) {
        Nv[i] = Nmin + i;
    }

    int numthreads = stoi(argv[9]);

    int L = 50;
    int nmax = 7;

    int nsweeps = 20;
    Real errgoal = -1;
    bool quiet = true;

    Sweeps sweeps(nsweeps);
    sweeps.minm() = 20;
    sweeps.maxm() = 10,20,100,100,200,200,300,300,400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 4;
    sweeps.noise() = 1e-3;//,1e-3,1-3,1e-3,1e-3,1e-4,1e-4,1e-4,1e-5,1e-5,1e-6,1e-7,1e-8,0;//1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,0.0;

    vector<int> minm(nsweeps);
    vector<int> maxm(nsweeps);
    vector<Real> cutoff(nsweeps);
    vector<int> niter(nsweeps);
    vector<Real> noise(nsweeps);
    for(int sw = 1; sw <= nsweeps; ++sw) {
        minm[sw-1] = sweeps.minm(sw);
        maxm[sw-1] = sweeps.maxm(sw);
        cutoff[sw-1] = sweeps.cutoff(sw);
        niter[sw-1] = sweeps.niter(sw);
        noise[sw-1] = sweeps.noise(sw);
    }

#ifdef MACOSX
    string resdir = "/Users/Abuenameh/Documents/Simulation Results/BH-ITensor-DMRG";
#endif
#ifdef AMAZON_EC2
    string resdir = "/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/BH-ITensor-DMRG";
#endif
    string resfile = format("%s/res.%d.txt", resdir, resi);
    ofstream os(resfile);

    printMath(os, "seed", resi, seed);
    printMath(os, "ximax", resi, ximax);
    printMath(os, "L", resi, L);
    printMath(os, "nmax", resi, nmax);
    printMath(os, "errgoal", resi, errgoal);
    printMath(os, "nsweep", resi, nsweeps);
    printMath(os, "maxm", resi, maxm);
    printMath(os, "cutoff", resi, cutoff);
    printMath(os, "niter", resi, niter);
    printMath(os, "noise", resi, noise);
    printMath(os, "ts", resi, xv);
    printMath(os, "Ns", resi, Nv);

    BoseHubbardSiteSet sites(L, nmax);

#ifdef MACOSX
    string groundstate = "/Users/Abuenameh/Projects/ITensorDMRG/GroundState/Release/groundstate";
#endif
#ifdef AMAZON_EC2
    string groundstate = "/home/ubuntu/ITensorDMRG/GroundState/Release/GroundState";
#endif
    ProcessPool pool(numthreads, groundstate, [&] (ostream& os, istream& is, istream& abortis, bool& abort) {
        write(os, sites);
        write(os, nsweeps);
        write(os, minm);
        write(os, maxm);
        write(os, niter);
        write(os, cutoff);
        write(os, noise);
        write(os, errgoal);
        write(os, quiet);
    });

    concurrent_queue<Results> resq;

    vector<Real> Us(L, 1);
    
    vector<Real> mus(L, 0);
    std::mt19937 gen;
    std::uniform_real_distribution<> dist(-ximax, ximax);
    auto randmu = bind(dist, gen);
    gen.seed(seed);
    for(int i = 0; i < L; i++) {
        mus[i] = randmu();
    }

    for(int ix = 0; ix < nx; ++ix) {
        for(int iN = 0; iN < nN; ++iN) {
            vector<Real> xs(L, xv[ix]);
            pool.enqueue([&](ostream& os, istream& is, istream& abortis, bool& abort, concurrent_queue<Results>* resq, int ix, int iN, vector<Real>& xs, vector<Real>& Us, vector<Real>& mus, int N) {
                write(os, xs);
                write(os, Us);
                write(os, mus);
                write(os, N);
                
                read(abortis, abort);
                if(abort) {
                    return;
                }

                Results res;

                res.ix = ix;
                res.iN = iN;
                res.xs = xs;
                res.Us = Us;
                res.mus = mus;

                read(is, res.E0);
                read(is, res.Ei);
                read(is, res.n);
                read(is, res.n2);
                read(is, res.C);
                read(is, res.runtime);

                resq->push(res);
                
            }, &resq, ix, iN, xs, Us, mus, Nv[iN]);
        }
    }

#ifdef MACOSX
    string python = "/Library/Frameworks/Python.framework/Versions/2.7/bin/python";
    string script = "/Users/Abuenameh/PycharmProjects/DMRG/ZMQProgressDialog.py";
#endif
#ifdef AMAZON_EC2
    string python = "/usr/bin/python";
    string script = "/home/ubuntu/PycharmProjects/BH-DMRG/ZMQProgressDialog.py";
#endif
    vector<string> args;
    args.push_back(python);
    args.push_back(script);

    child c = execute(set_args(args), inherit_env());

    zmq::context_t context(1);
    zmq::socket_t socket(context, ZMQ_PUSH);

    socket.connect("tcp://localhost:5556");
    string len = to_string(nx*nN);
    zmq::message_t lenmessage(len.length());
    memcpy ((void *) lenmessage.data(), len.c_str(), len.length());
    socket.send(lenmessage);

    interrupted = false;
    signal(SIGINT, sigint);

    int count = 0;

    multi_array<vector<Real>, 2> xres(extents[nx][nN]);
    multi_array<vector<Real>, 2> Ures(extents[nx][nN]);
    multi_array<vector<Real>, 2> mures(extents[nx][nN]);
    multi_array<Real, 2> E0res(extents[nx][nN]);
    multi_array<vector<Real>, 2> Eires(extents[nx][nN]);
    multi_array<vector<Real>, 2> nres(extents[nx][nN]);
    multi_array<vector<Real>, 2> n2res(extents[nx][nN]);
    multi_array<vector<vector<Real> >, 2> Cres(extents[nx][nN]);
    multi_array<string, 2> runtimei(extents[nx][nN]);

    for(int ix = 0; ix < nx; ++ix) {
        for(int iN = 0; iN < nN; ++iN) {
            xres[ix][iN] = vector<Real>(L, NAN);
            Ures[ix][iN] = vector<Real>(L, NAN);
            mures[ix][iN] = vector<Real>(L, NAN);
            E0res[ix][iN] = NAN;
            Eires[ix][iN] = vector<Real>(1, NAN);
            nres[ix][iN].assign(L, NAN);
            n2res[ix][iN].assign(L, NAN);
            Cres[ix][iN].assign(L, vector<Real>(L, NAN));
            runtimei[ix][iN] = "unfinished";
        }
    }

    Results res;
    while(!interrupted && count++ < nx*nN) {
        resq.wait_and_pop(res);
        int ix = res.ix;
        int iN = res.iN;
        xres[ix][iN] = res.xs;
        Ures[ix][iN] = res.Us;
        mures[ix][iN] = res.mus;
        E0res[ix][iN] = res.E0;
        Eires[ix][iN] = res.Ei;
        
        nres[ix][iN] = res.n;
        n2res[ix][iN] = res.n2;
        Cres[ix][iN] = res.C;

        runtimei[ix][iN] = seconds_to_string(res.runtime);

        zmq::message_t message(sizeof(int));
        ((int*)message.data())[0] = 1;
        socket.send(message);
    }

    pool.interrupt();
    terminate(c);

    printMath(os, "tres", resi, xres);
    printMath(os, "Ures", resi, Ures);
    printMath(os, "mures", resi, mures);
    printMath(os, "E0res", resi, E0res);
    printMath(os, "Eires", resi, Eires);
    printMath(os, "nres", resi, nres);
    printMath(os, "n2res", resi, n2res);
    printMath(os, "Cres", resi, Cres);
    printMath(os, "runtimei", resi, runtimei);

    time_point<system_clock> end = system_clock::now();
    string runtime = seconds_to_string(duration_cast<seconds>(end - start).count());
    printMath(os, "runtime", resi, runtime);

    exit(0);

    return 0;
}
