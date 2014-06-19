#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>

#include <boost/multi_array.hpp>
#include <boost/process.hpp>

#include <zmq.hpp>

#undef small
#include <core.h>
#include <hambuilder.h>

#include "ipc.h"
#include "concurrent_queue.h"
#include "ThreadPool.h"
#include "mathematica.h"
#include "ProcessPool.h"

#include "BoseHubbardSiteSet.h"
#include "BoseHubbardHamiltonian.h"
#include "BoseHubbardObserver.h"

#include <sites/hubbard.h>
#include <hams/HubbardChain.h>

using std::ref;
using std::stoi;
using std::stod;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
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

using namespace zmq;

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
    //vector<vector<Real> > C;
    vector<Real> C;
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

Results* interruptRes;
concurrent_queue<Results>* interruptResq;

#ifdef FST

BOOL WINAPI ConsoleHandler(DWORD type)
{
    if(type == CTRL_C_EVENT) {
        interrupted = true;
        interruptResq->push(*interruptRes);
        return TRUE;
    }
    return FALSE;
}

#else

void sigint(int)
{
    interrupted = true;
    interruptResq->push(*interruptRes);
    cout << "Interrupted" << endl;
}

#endif

int main(int argc, char **argv)
{
    if(argc < 10) {
        cerr << "Insufficient number of arguments" << endl;
        return 1;
    }

    time_point<system_clock> start = system_clock::now();
    
    ProcessPool::queue_idx = start.time_since_epoch().count();

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

    int L = 5;
    int nmax = 7;

    int nsweeps = 40;
    Real errgoal = -1;
    bool quiet = true;

    Sweeps sweeps(nsweeps);
    sweeps.minm() = 20;
    sweeps.maxm() = 400;//10,20,100,100,200,200,300,300,400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 4;
    sweeps.noise() = 1,1,1,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,0;//1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,0;//1e-3;//,1e-3,1-3,1e-3,1e-3,1e-4,1e-4,1e-4,1e-5,1e-5,1e-6,1e-7,1e-8,0;//1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,0.0;

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
//    string resdir = "/Users/Abuenameh/Documents/Simulation Results/BH-ITensor-DMRG";
    string resdir = "/Users/Abuenameh/Dropbox/Simulation Results/ITensorDMRG";
#endif
#ifdef AMAZON_EC2
    string resdir = "/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/BH-ITensor-DMRG";
#endif
#ifdef FST
    string resdir = "C:/Users/abuenameh/Dropbox/Server/ITensorDMRG";
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
    
    vector<IQMPO> Cds;
    vector<MPO> Cdstmp;

    string setupfile = format("%s/setup.%d.%d.dat", resdir, L, nmax);
    ifstream setupis(setupfile+"a", ios::binary);
    if(setupis.good()) {
        read(setupis, sites);
        Cdstmp.resize(L-1, MPO(sites));
        read(setupis, Cdstmp);
    } else {
        for(int d = 1; d < L; d++) {
            MPO Cd = HamBuilder<ITensor>(sites, "Bdag", 1, "B", 1+d);
            for(int i = 2; i <= L-d; i++) {
                Cd.plusEq(HamBuilder<ITensor>(sites, "Bdag", i, "B", i+d));
            }
            Cd *= 1./(L-d);
            Cdstmp.push_back(Cd);
            //Cds.push_back(Cd.toIQMPO());
        }
        ofstream setupos(setupfile, ios::binary);
        write(setupos, sites);
        write(setupos, Cdstmp);
    }
    for(int i = 0; i < L-1; i++) {
        Cds.push_back(Cdstmp[i].toIQMPO());
    }
    
    /*context_t gscontext(numthreads);
    vector<context_t> incontexts;
    vector<socket_t> outsockets, insockets;
    vector<connection_monitor> monitors;
    vector<thread> monitorthreads;
    string monitoraddr = "inproc://monitor.";
//    vector<int> outports, inports;
//    int port = 5555;
    for(int i = 0; i < numthreads; i++) {
        outsockets.emplace_back(gscontext, ZMQ_PUSH);
//        outports.push_back(port++);
        //outsockets.back().connect(("tcp://localhost:" + to_string(outports.back())).c_str());

        incontexts.emplace_back();
        insockets.emplace_back(incontexts.back(), ZMQ_PULL);
        monitors.emplace_back(incontexts.back(), insockets.back());
        monitorthreads.emplace_back([&] () {
            monitors.back().monitor(insockets.back(), (monitoraddr + to_string(reinterpret_cast<int>(&insockets.back()))).c_str(), ZMQ_EVENT_DISCONNECTED);
        });
//        inports.push_back(port++);
        //insockets.back().bind(("tcp://localhost:" + to_string(inports.back())).c_str());
    }*/
    
    
#ifdef MACOSX
//    string groundstate = "/Users/Abuenameh/Projects/ITensorDMRG/GroundState/Release/groundstate";
    string groundstate = "/Users/Abuenameh/NetBeansProjects/DMRGGroundState/dist/Release/CLang-MacOSX/dmrggroundstate";
#endif
#ifdef AMAZON_EC2
    string groundstate = "/home/ubuntu/ITensorDMRG/GroundState/Release/GroundState";
#endif
#ifdef FST
    string groundstate = "C:/Users/abuenameh/Documents/NetBeansProjects/DMRGGroundState/dist/Release/MinGW_TDM-Windows/dmrggroundstate.exe";
#endif
    ProcessPool pool(numthreads, groundstate, /*outsockets, incontexts, insockets,*/ [&] (message_queue& mq,/*nnxx::socket& os, nnxx::socket& is,*/ bool& abort) {
        cout << "Writing sites" << endl;
        write(mq, sites);
        cout << "Writing Cds" << endl << flush;
        write(mq, Cds);
        cout << "Wrote Cds" << endl << flush;
        write(mq, nsweeps);
        cout << "Wrote nsweeps" << endl << flush;
        write(mq, minm);
        write(mq, maxm);
        write(mq, niter);
        write(mq, cutoff);
        write(mq, noise);
        cout << "Wrote sweeps" << endl << flush;
        write(mq, errgoal);
        cout << "Wrote errgoal" << endl << flush;
        write(mq, quiet);
        cout << "Wrote quiet" << endl << flush;
        cout << "About to wait" << endl << flush;
        int wait = 0;
        read(mq, wait);
    });
    
    concurrent_queue<Results> resq;
    
//#ifndef FST
    interruptResq = &resq;
    interruptRes = new Results;
//#endif
    
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
            pool.enqueue([&](message_queue& mq,/*nnxx::socket& os, nnxx::socket& is,*/ bool& abort, concurrent_queue<Results>* resq, int ix, int iN, vector<Real>& xs, vector<Real>& Us, vector<Real>& mus, int N) {
                cout << "Writing parameters " << (int)(&mq) << endl << flush;
//                for(;;) {}
                try{
                write(mq, xs);
                cout << "Wrote xs " << (int)(&mq) << endl << flush;
                write(mq, Us);
                cout << "Wrote Us " << (int)(&mq) << endl << flush;
                write(mq, mus);
                cout << "Wrote mus" << endl << flush;
                write(mq, N);
                cout << "Wrote N" << endl << flush;
                }catch(std::exception& e) {cout << "Error: " << e.what() << endl << flush; abort=true; return; }

#ifndef FST
                read(abortis, abort);
                if(abort) {
                    return;
                }
#endif
                
                Results res;

                res.ix = ix;
                res.iN = iN;
                res.xs = xs;
                res.Us = Us;
                res.mus = mus;

                try{
                    cout << "About to read E0" << endl << flush;
                read(mq, res.E0);
                read(mq, res.Ei);
                read(mq, res.n);
                read(mq, res.n2);
                read(mq, res.C);
                read(mq, res.runtime);

                resq->push(res);
                }
                catch(...) {
                    cout << "-----Aborting-----" << endl << endl;
                    abort = true;
                    return;
                }
                abort = false;

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
#ifdef FST
    string python = "C:/Python27/python.exe";
    string script = "C:/Users/abuenameh/PycharmProjects/BH-DMRG/ZMQProgressDialog.py";
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
#ifdef FST
    SetConsoleCtrlHandler((PHANDLER_ROUTINE)ConsoleHandler, TRUE);
#else
    signal(SIGINT, sigint);
#endif
    
    int count = 0;

    multi_array<vector<Real>, 2> xres(extents[nx][nN]);
    multi_array<vector<Real>, 2> Ures(extents[nx][nN]);
    multi_array<vector<Real>, 2> mures(extents[nx][nN]);
    multi_array<Real, 2> E0res(extents[nx][nN]);
    multi_array<vector<Real>, 2> Eires(extents[nx][nN]);
    multi_array<vector<Real>, 2> nres(extents[nx][nN]);
    multi_array<vector<Real>, 2> n2res(extents[nx][nN]);
    //multi_array<vector<vector<Real> >, 2> Cres(extents[nx][nN]);
    multi_array<vector<Real>, 2> Cres(extents[nx][nN]);
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
            //Cres[ix][iN].assign(L, vector<Real>(L, NAN));
            Cres[ix][iN].assign(L, NAN);
            runtimei[ix][iN] = "unfinished";
        }
    }

    Results res;
    while(!interrupted && count++ < nx*nN) {
        resq.wait_and_pop(res);
        if(interrupted) {
            break;
        }
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
    try {
    terminate(c);
    }
    catch(...) {}

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
