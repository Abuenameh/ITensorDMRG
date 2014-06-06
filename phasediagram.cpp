#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <core.h>
#include <hambuilder.h>

#include <boost/multi_array.hpp>
#include <boost/process.hpp>

#include <zmq.hpp>

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
    }
    else {
        return format("%d:%02d:%02d", h, m, s);
    }
}



        /*for(int j = 1; j <= L; ++j) {
            psi.position(j);
            n(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("N",j)*psi.A(j));
            n2(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("N2",j)*psi.A(j));
            for(int k = 1; k <= L; ++k) {
                C(j, k) = psiHphi(psi,COp[j-1][k-1],psi);
            }
        }


        time_point<system_clock> end = system_clock::now();
        //seconds runtime = duration_cast<seconds>(end - start);
        //res.runtime = runtime.count();
        runtime = seconds_to_string(duration_cast<seconds>(end - start).count());*/


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

    const Real tmin = stod(argv[4]);
    const Real tmax = stod(argv[5]);
    const int nt = stoi(argv[6]);
    vector<Real> tv = linspace(tmin, tmax, nt);

    const int Nmin = stoi(argv[7]);
    const int Nmax = stoi(argv[8]);
    const int nN = Nmax - Nmin + 1;
    vector<int> Nv(nN);
    for(int i = 0; i < nN; ++i) {
        Nv[i] = Nmin + i;
    }

    int numthreads = stoi(argv[9]);

    int L = 20;
    int nmax = 7;

    int nsweeps = 5;
    Real errgoal = -1;
    bool quiet = true;

    Sweeps sweeps(nsweeps);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    
    vector<int> maxm(nsweeps);
    vector<Real> cutoff(nsweeps);
    vector<int> niter(nsweeps);
    vector<Real> noise(nsweeps);
    for(int sw = 0; sw < nsweeps; ++sw) {
        maxm[sw] = sweeps.maxm(sw);
        cutoff[sw] = sweeps.cutoff(sw);
        niter[sw] = sweeps.niter(sw);
        noise[sw] = sweeps.noise(sw);
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
    printMath(os, "ts", resi, tv);
    printMath(os, "Ns", resi, Nv);
    
    BoseHubbardSiteSet sites;

    vector<vector<IQMPO> > COp;

    try {

        sites = BoseHubbardSiteSet(L, nmax);

        vector<IQIndex> links(L+1);
        for(int l = 0; l <= L; ++l) {
            vector<IndexQN> indices;
            for(int n = 0; n <= nmax; n++) {
                indices.push_back(IndexQN(Index(nameint("n = ", n) + nameint(" for site ", l),1), QN(0,n,n%2)));
            }
            links.at(l) = IQIndex(nameint("BoseHubbard site=",l),indices);
        }

        /*for(int i = 1; i <= L; ++i) {
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
        }*/

    } catch(ITError e) {
        cerr << "ITensor error: " << e.what() << endl;
        return 1;
    } catch(...) {
        cerr << "Unknown error" << endl;
        return 1;
    }
    
    ProcessPool pool(numthreads, "/Users/Abuenameh/Projects/ITensorDMRG/GroundState/Release/groundstate", [&sites] (socket_t& out, socket_t& in, ostream& os) {
        cout << "Init called" << endl;
        stringstream ss(stringstream::in | stringstream::out | stringstream::binary);
        sites.write(ss);
        //out.send(ss.str().data(), ss.str().size());
        sites.write(os);
        cout << "Init finished" << endl;
    });

    /*ThreadPool pool(numthreads);
    concurrent_queue<Results> resq;
    concurrent_queue<int> q;

    vector<Real> Us(L, 1);
    vector<Real> mus(L, 0);

    multi_array<vector<Real>, 2> tres(extents[nt][nN]);
    multi_array<vector<Real>, 2> Ures(extents[nt][nN]);
    multi_array<vector<Real>, 2> mures(extents[nt][nN]);
    multi_array<Real, 2> E0res(extents[nt][nN]);
    multi_array<vector<Real>, 2> Eires(extents[nt][nN]);
    multi_array<Vector, 2> nres(extents[nt][nN]);
    multi_array<Vector, 2> n2res(extents[nt][nN]);
    multi_array<Matrix, 2> Cres(extents[nt][nN]);
    multi_array<string, 2> runtimei(extents[nt][nN]);

    for(int it = 0; it < nt; ++it) {
        for(int iN = 0; iN < nN; ++iN) {
            tres[it][iN] = vector<Real>(L, NAN);
            Ures[it][iN] = vector<Real>(L, NAN);
            mures[it][iN] = vector<Real>(L, NAN);
            E0res[it][iN] = NAN;
            Eires[it][iN] = vector<Real>(1, NAN);
            nres[it][iN] = Vector(L);
            nres[it][iN] = NAN;
            n2res[it][iN] = Vector(L);
            n2res[it][iN] = NAN;
            Cres[it][iN] = Matrix(L, L);
            for(int i = 1; i <= L; ++i) {
                Cres[it][iN].Row(i) = NAN;
            }
            runtimei[it][iN] = "unfinished";
        }
    }
    
    for(int it = 0; it < nt; ++it) {
        for(int iN = 0; iN < nN; ++iN) {
            vector<Real> ts(L, tv[it]);
            tres[it][iN] = ts;
            Ures[it][iN] = Us;
            mures[it][iN] = mus;
            pool.enqueue(bind(groundstate, ref(q), ref(sweeps), errgoal, quiet, ref(sites), ref(COp), it, iN, ts, Us, mus, Nv[iN], ref(E0res[it][iN]), ref(Eires[it][iN]), ref(nres[it][iN]), ref(n2res[it][iN]), ref(Cres[it][iN]), ref(runtimei[it][iN])));
        }
    }*/
    
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
    string len = to_string(nt*nN);
    zmq::message_t lenmessage(len.length());
    memcpy ((void *) lenmessage.data(), len.c_str(), len.length());
    socket.send(lenmessage);

    interrupted = false;
    signal(SIGINT, sigint);

    int count = 0;

    //Results res;
    /*multi_array<vector<Real>, 2> tres(extents[nt][nN]);
    multi_array<vector<Real>, 2> Ures(extents[nt][nN]);
    multi_array<vector<Real>, 2> mures(extents[nt][nN]);
    multi_array<Real, 2> E0res(extents[nt][nN]);
    multi_array<vector<Real>, 2> Eires(extents[nt][nN]);
    multi_array<Vector, 2> nres(extents[nt][nN]);
    multi_array<Vector, 2> n2res(extents[nt][nN]);
    multi_array<Matrix, 2> Cres(extents[nt][nN]);
    multi_array<string, 2> runtimei(extents[nt][nN]);*/
    
    /*for(int it = 0; it < nt; ++it) {
        for(int iN = 0; iN < nN; ++iN) {
            tres[it][iN] = vector<Real>(L, NAN);
            Ures[it][iN] = vector<Real>(L, NAN);
            mures[it][iN] = vector<Real>(L, NAN);
            E0res[it][iN] = NAN;
            Eires[it][iN] = vector<Real>(1, NAN);
            nres[it][iN] = Vector(L);
            nres[it][iN] = NAN;
            n2res[it][iN] = Vector(L);
            n2res[it][iN] = NAN;
            Cres[it][iN] = Matrix(L, L);
            for(int i = 1; i <= L; ++i) {
                Cres[it][iN].Row(i) = NAN;
            }
            runtimei[it][iN] = "unfinished";
        }
    }*/
    
    int qi;
    while(!interrupted && count++ < nt*nN) {
        //q.wait_and_pop(qi);
        //resq.wait_and_pop(res);
        /*int it = res.it;
        int iN = res.iN;
        tres[it][iN] = res.ts;
        Ures[it][iN] = res.Us;
        mures[it][iN] = res.mus;
        E0res[it][iN] = res.E0;
        Eires[it][iN] = res.Ei;
        nres[it][iN] = res.n;
        n2res[it][iN] = res.n2;
        Cres[it][iN] = res.C;
        runtimei[it][iN] = seconds_to_string(res.runtime);*/
        zmq::message_t message(sizeof(int));
        ((int*)message.data())[0] = 1;
        socket.send(message);
    }
    
    //pool.interrupt();
    terminate(c);

    /*printMath(os, "tres", resi, tres);
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
    printMath(os, "runtime", resi, runtime);*/
	
    exit(0);
	
    return 0;
}
