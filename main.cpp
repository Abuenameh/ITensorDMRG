#include <iostream>
#include <fstream>
#include <chrono>
#include <core.h>
#include <hambuilder.h>

#include "BoseHubbardSiteSet.h"
#include "BoseHubbardHamiltonian.h"
#include "BoseHubbardObserver.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::chrono::time_point;
using std::chrono::duration;
using std::chrono::system_clock;

using namespace itensor;

int main(int argc, char **argv)
{
    if(argc != 5) {
        cerr << "Missing parameter or output file name" << endl;
        return 1;
    }

    string command(argv[1]);
    if(command != "setup" && command != "run") {
        cerr << "Unknown command" << endl;
        return 1;
    }

    if(command == "setup") {
        ofstream setup(argv[2]);

        const int L = std::stoi(argv[3]);
        const int nmax = std::stoi(argv[4]);
        BoseHubbardSiteSet sites(L, nmax);
        sites.write(setup);

        std::vector<IQIndex> links(L+1);
        for(int l = 0; l <= L; ++l) {
            std::vector<IndexQN> indices;
            for(int n = 0; n <= nmax; n++) {
                indices.push_back(IndexQN(Index(nameint("n = ", n) + nameint(" for site ", l),1), QN(0,n,n%2)));
            }
            links.at(l) = IQIndex(nameint("BoseHubbard site=",l),indices);
        }

        for(int i = 1; i <= L; ++i) {
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
                Cij.write(setup);
            }
        }
    } else if(command == "run") {
        ifstream setup(argv[2]);


        time_point<system_clock> start = system_clock::now();

        string parmfilename(argv[3]);
        InputFile parmfile(parmfilename);
        InputGroup parms(parmfile,"parameters");

        const int L = parms.getInt("L");
        const int nmax = parms.getInt("nmax");
        const int nsweeps = parms.getInt("nsweeps");
        const Real errgoal = parms.getReal("errgoal");
        const int N = parms.getInt("N");
        const bool quiet = parms.getYesNo("quiet",true);

        InputGroup table(parms,"sweeps");
        Sweeps sweeps(nsweeps,table);

        string t = parms.getString("t", "0.01");
        string U = parms.getString("U", "1");
        string mu = parms.getString("mu", "0");

        BoseHubbardSiteSet sites;
        sites.read(setup);

        BoseHubbardHamiltonian BH = BoseHubbardHamiltonian(sites, Opt("t",t)&Opt("U",U)&Opt("mu",mu));
        IQMPO H = BH;

        vector<Real> ts = BH.t();
        vector<Real> Us = BH.U();
        vector<Real> mus = BH.mu();

        InitState initState(sites);
        int n0 = N / L;
        for(int i = 1; i <= L; ++i) {
            if(i <= N % L)
                initState.set(i,std::to_string(n0+1));
            else
                initState.set(i,std::to_string(n0));
        }
        IQMPS psi(initState);

        vector<vector<IQMPO> > COp;

        for(int i = 1; i <= L; ++i) {
            std::vector<IQMPO> Ci;
            for(int j = 1; j <= L; ++j) {
                IQMPO Cij(sites);
                Cij.read(setup);
                Ci.push_back(Cij);
            }
            COp.push_back(Ci);
        }

        /*std::vector<IQIndex> links(L+1);
        for(int l = 0; l <= L; ++l) {
            std::vector<IndexQN> indices;
            for(int n = 0; n <= nmax; n++) {
                indices.push_back(IndexQN(Index(nameint("n = ", n) + nameint(" for site ", l),1), QN(0,n,n%2)));
            }
            links.at(l) = IQIndex(nameint("BoseHubbard site=",l),indices);
        }

        vector<vector<IQMPO> > COp;

        for(int i = 1; i <= L; ++i) {
            std::vector<IQMPO> Ci;
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

        std::cout << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << std::endl;

        std::cout << sweeps;

        BoseHubbardObserver<IQTensor> observer(psi,Opt("EnergyErrgoal",errgoal));

        Real E0 = dmrg(psi,H,sweeps,observer,Opt("Quiet",quiet));

        std::cout << format("\nGround State Energy = %.10f",E0) << std::endl;

        std::vector<Real> Ei = observer.getEnergies();

        Real Ntot = 0;

        Vector n(L), n2(L);
        Matrix C(L, L);
        for(int j = 1; j <= L; ++j) {
            psi.position(j);
            n(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("N",j)*psi.A(j));
            n2(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("N2",j)*psi.A(j));
            Ntot += n(j);
            for(int k = 1; k <= L; ++k) {
                C(j, k) = psiHphi(psi,COp[j-1][k-1],psi);
            }
        }

        string outputfilename(argv[3]);
        ofstream outputfile(outputfilename);

        outputfile << "t ";
        for(unsigned int j = 0; j < ts.size(); ++j) {
            outputfile << format("%.20e ", ts[j]);
        }
        outputfile << endl;
        outputfile << "U ";
        for(unsigned int j = 0; j < Us.size(); ++j) {
            outputfile << format("%.20e ", Us[j]);
        }
        outputfile << endl;
        outputfile << "mu ";
        for(unsigned int j = 0; j < mus.size(); ++j) {
            outputfile << format("%.20e ", mus[j]);
        }
        outputfile << endl;
        outputfile << "Ei ";
        for(unsigned int j = 0; j < Ei.size(); ++j) {
            outputfile << format("%.20e ", Ei[j]);
        }
        outputfile << endl;
        outputfile << "E0 " << format("%.20e", E0) << endl;
        outputfile << "n ";
        for(int j = 1; j <= L; ++j) {
            outputfile << format("%.10e ", n(j));
        }
        outputfile << endl;
        outputfile << "n2 ";
        for(int j = 1; j <= L; ++j) {
            outputfile << format("%.10e ", n2(j));
        }
        outputfile << endl;
        outputfile << "C ";
        for(int j = 1; j <= L; ++j) {
            for(int k = 1; k <= L; ++k) {
                outputfile << format("%.10e ", C(j,k));
            }
        }
        outputfile << endl;

        time_point<system_clock> end = system_clock::now();
        duration<double> runtime = end - start;
        outputfile << "runtime " << runtime.count() << endl;

        /*cout << "Density:" << endl;
        for(int j = 1; j <= L; ++j)
            cout << format("%d %.10f\n", j, n(j));
        cout << endl;

        cout << "Density^2:" << endl;
        for(int j = 1; j <= L; ++j)
            cout << format("%d %.10f\n", j, n2(j));
        cout << endl;

        cout << endl << format("Ntot = %.10f\n", Ntot) << endl;

        cout << "Correlation:" << endl;
        for(int j = 1; j <= L; ++j) {
            for(int k = 1; k <= L; ++k)
                cout << format("%.10f ", C(j,k));
            cout << endl;
        }*/
    }

    return 0;
}
