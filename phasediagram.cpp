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
    //IQMPS psi0;
    string psi0;
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

    int L = 20;
    int nmax = 7;

    int nsweeps = 5;
    Real errgoal = -1;
    bool quiet = true;

    Sweeps sweeps(nsweeps);
    sweeps.minm() = 20;
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;

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

    BoseHubbardSiteSet sites;

    vector<vector<IQMPO> > COp;
    vector<IQMPO> BOp;
    
    vector<IQMPO> bs;

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
            IQMPO Bi = HamBuilder<IQTensor>(sites, "B", i);
            Bi.position(1);
            /*for(int n = 1; n <= L; ++n) {
                IQTensor& W = Bi.Anc(n);
                IQIndex row = dag(links[n-1]), col = links[n];

                W = IQTensor(dag(sites.si(n)),sites.siP(n),row,col);

            W += sites.op("Id",n) * row(1) * col(1);
            W += sites.op("Id",n) * row(2) * col(2);
            W += sites.op("Id",n) * row(3) * col(3);
            W += sites.op("Id",n) * row(4) * col(4);

                if(n == i) {
                    W += sites.op("B",n) * row(2) * col(4);
                } else {
                    for(int k = 1; k <= 4; ++k) {
                        W += sites.op("Id",n) * row(k) * col(k);
                    }
                }
            }
            //Bi.Anc(1) *= IQTensor(links.at(0)(1));
            //Bi.Anc(L) *= IQTensor(dag(links.at(L)(4)));*/
            BOp.push_back(Bi);
            COp.push_back(Ci);
            bs.push_back(Bi);
        }

    } catch(ITError e) {
        cerr << "ITensor error: " << e.what() << endl;
        return 1;
    } catch(...) {
        cerr << "Unknown error" << endl;
        return 1;
    }

#ifdef MACOSX
    string groundstate = "/Users/Abuenameh/Projects/ITensorDMRG/GroundState/Release/groundstate";
#endif
#ifdef AMAZON_EC2
    string groundstate = "/home/ubuntu/ITensorDMRG/GroundState/Release/GroundState";
#endif
    ProcessPool pool(numthreads, groundstate, [&] (ostream& os, istream& is) {
        write(os, sites);
        int test = 42;
        //write(os, test);
        IQMPS m(sites);
        //write(os, m);
        //cout << m << endl << endl << endl;
        IQMPO b = bs[0];
        write(os, b);
        os.flush();
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

    for(int ix = 0; ix < nx; ++ix) {
        for(int iN = 0; iN < nN; ++iN) {
            vector<Real> xs(L, xv[ix]);
            pool.enqueue([&](ostream& os, istream& is, concurrent_queue<Results>* resq, int ix, int iN, vector<Real>& xs, vector<Real>& Us, vector<Real>& mus, int N) {
                write(os, xs);
                write(os, Us);
                write(os, mus);
                write(os, N);

                Results res;

                res.ix = ix;
                res.iN = iN;
                res.xs = xs;
                res.Us = Us;
                res.mus = mus;

                string psi0;
                int len;
                read(is, len);
                vector<char> buf(len);
                is.read(buf.data(), len);
                psi0.append(buf.data(), len);
                res.psi0 = psi0;

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

    zmq::context_t context;
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
    multi_array<Vector, 2> nres(extents[nx][nN]);
    multi_array<Vector, 2> n2res(extents[nx][nN]);
    multi_array<Matrix, 2> Cres(extents[nx][nN]);
    multi_array<string, 2> runtimei(extents[nx][nN]);

    for(int ix = 0; ix < nx; ++ix) {
        for(int iN = 0; iN < nN; ++iN) {
            xres[ix][iN] = vector<Real>(L, NAN);
            Ures[ix][iN] = vector<Real>(L, NAN);
            mures[ix][iN] = vector<Real>(L, NAN);
            E0res[ix][iN] = NAN;
            Eires[ix][iN] = vector<Real>(1, NAN);
            nres[ix][iN] = Vector(L);
            nres[ix][iN] = NAN;
            n2res[ix][iN] = Vector(L);
            n2res[ix][iN] = NAN;
            Cres[ix][iN] = Matrix(L, L);
            for(int i = 1; i <= L; ++i) {
                Cres[ix][iN].Row(i) = NAN;
            }
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

        stringstream ss(res.psi0, std::ios_base::in);
        IQMPS psi0(sites);
        psi0.read(ss);
        bool old = false;
        bool exact = true;
        vector<IQMPS> Bpsi;
        psi0.position(1);
        time_point<system_clock> startapply = system_clock::now();
        if(!old) {
            for(int j = 1; j <= L; ++j) {
                IQMPS psi;
                if(exact)
                    exactApplyMPO(psi0, BOp[j-1], psi);
                else
                    zipUpApplyMPO(psi0, BOp[j-1], psi);
                Bpsi.push_back(psi);
            }
        }
        time_point<system_clock> endapply = system_clock::now();
        //seconds runtime = duration_cast<seconds>(end - start);
        //res.runtime = runtime.count();
        //cout <<  seconds_to_string(duration_cast<seconds>(endapply - startapply).count()) << endl;
        for(int j = 1; j <= L; ++j) {
            psi0.position(j);
            nres[ix][iN](j) = Dot(conj(primed(psi0.A(j),Site)),sites.op("N",j)*psi0.A(j));
            n2res[ix][iN](j) = Dot(conj(primed(psi0.A(j),Site)),sites.op("N2",j)*psi0.A(j));
            for(int k = 1; k <= L; ++k) {
                if(old)
                    Cres[ix][iN](j, k) = psiHphi(psi0,COp[j-1][k-1],psi0);
                else {
                    /*IQMPS psiL, psiR;
                    if(exact) {
                        exactApplyMPO(psi0, BOp[j-1], psiL);
                        exactApplyMPO(psi0, BOp[k-1], psiR);

                    } else {
                        zipUpApplyMPO(psi0, BOp[j-1], psiL);
                        zipUpApplyMPO(psi0, BOp[k-1], psiR);

                    }
                    Cres[ix][iN](j, k) = psiphi(psiL, psiR);*/
                    Cres[ix][iN](j, k) = psiphi(Bpsi[j-1], Bpsi[k-1]);

                }
            }
        }
        /*HamBuilder<IQTensor> hb(sites, "B", 1);
        IQMPO qwe = hb;//BOp[0];// = sites.op("B", 1);
        //cout << qwe << endl;
        qwe.position(1);
        psi0.position(1);
        IQMPS wer;// = qwe*psi0;
        //zipUpApplyMPO(psi0,qwe,wer);
        exactApplyMPO(psi0,qwe,wer);
        //cout << wer << endl;

        IQMPO asd = HamBuilder<IQTensor>(sites, "B", 2);//BOp[1];//sites.op("B", 2);
        asd.position(1);
        psi0.position(1);
        IQMPS sdf;
        //zipUpApplyMPO(psi0,asd,sdf);
        exactApplyMPO(psi0,asd,sdf);

        Real zxc = psiphi(wer,wer);
        cout << "zxc = " << zxc << endl;
        Real xcv = psiphi(wer,sdf);
        cout << "xcv = " << xcv << endl;*/
        //IQMPS dfg = conj(primed(sdf,Site));
        //cout << dfg << endl;

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
