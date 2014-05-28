#include <iostream>
#include <boost/format.hpp>
#include <core.h>

#include "BoseHubbardModel.h"
#include "BoseHubbardHamiltonian.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using boost::format;

int main(int argc, char **argv)
{
    if(argc != 2) {
        cerr << "Missing parameter file name" << endl;
        return 1;
    }
    
    string parmfilename(argv[1]);
    InputFile parmfile(parmfilename);
    InputGroup parms(parmfile,"parameters");
    
    const int L = parms.getInt("L");
    const int nmax = parms.getInt("nmax");
    const int nsweeps = parms.getInt("nsweeps");
    const int N = parms.getInt("N");
    
    BoseHubbardModel model(L, nmax);
    IQMPO H = BoseHubbardHamiltonian(model);

    InitState initState(model);
    for(int i = 1; i <= L; ++i) {
        if(i <= N)
            initState.set(i,"1");
        else
            initState.set(i,"0");
    }
    IQMPS psi(initState);

    std::vector<Index> links(L+1);
    for(int l = 0; l <= L; ++l) links.at(l) = Index(nameint("hl",l),4);

    std::vector<std::vector<IQMPO> > C;

    for(int i = 1; i <= L; ++i) {
        std::vector<IQMPO> Ci;
        for(int j = 1; j <= L; ++j) {
            MPO Cij(model);
            for(int n = 1; n <= L; ++n) {
                ITensor& W = Cij.Anc(n);
                Index &row = links[n-1], &col = links[n];

                W = ITensor(model.si(n),model.siP(n),row,col);

                if(i != j && (n == i || n == j)) {
                    W += model.op("Bdag",n) * row(1) * col(2);
                    W += model.op("B",n) * row(2) * col(4);
                } else {
                    for(int k = 1; k <= 4; ++k) {
                        W += model.op("Id",n) * row(k) * col(k);
                    }
                }
            }
            Cij.Anc(1) *= ITensor(links.at(0)(1));
            Cij.Anc(L) *= ITensor(links.at(L)(4));
            Ci.push_back(Cij);
        }
        C.push_back(Ci);
    }

    std::cout << boost::format("Initial energy = %.5f") % psiHphi(psi,H,psi) << std::endl;

    Sweeps sweeps(10);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-6,1E-7,0.0;
    std::cout << sweeps;

    Real En = dmrg(psi,H,sweeps,Quiet());

    std::cout << boost::format("\nGround State Energy = %.10f")%En << std::endl;

    Real Ntot = 0;

    Vector nd(L), n2d(L);
    for(int j = 1; j <= L; ++j) {
        psi.position(j);
        nd(j) = Dot(conj(primed(psi.A(j),Site)),model.op("N",j)*psi.A(j));
        n2d(j) = Dot(conj(primed(psi.A(j),Site)),model.op("N2",j)*psi.A(j));
        Ntot += nd(j);
    }

    cout << "Density:" << endl;
    for(int j = 1; j <= L; ++j)
        cout << format("%d %.10f\n") % j % nd(j);
    cout << endl;

    cout << "Density^2:" << endl;
    for(int j = 1; j <= L; ++j)
        cout << format("%d %.10f\n") % j % n2d(j);
    cout << endl;

    cout << endl << format("Ntot = %.10f\n") % Ntot << endl;

    cout << "Correlation:" << endl;
    for(int j = 1; j <= L; ++j) {
        for(int k = 1; k <= L; ++k)
            cout << format("%.10f ") % psiHphi(psi,C[j-1][k-1],psi);
        cout << endl;
    }

    return 0;
}
