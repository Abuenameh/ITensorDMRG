#include <iostream>
#include <fstream>
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

using namespace itensor;

int main(int argc, char **argv)
{
    if(argc != 3) {
        cerr << "Missing parameter or output file name" << endl;
        return 1;
    }
    
    string parmfilename(argv[1]);
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
    string mu = parms.getString("mu", "0.0244067519637,0.107594683186,0.0513816880358,0.0224415914984,-0.0381726003305,0.0729470565333,-0.0312063943687,0.195886500391,0.231831380251,-0.0582792405871,0.145862519041,0.0144474598765");
    
    /*int L = 12;
    int nmax = 7;
    int nsweeps = 4;
    int N = 4;*/
    
    BoseHubbardSiteSet model(L, nmax);
    IQMPO H = BoseHubbardHamiltonian(model, Opt("t",t)&Opt("U",U)&Opt("mu",mu));

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

    vector<vector<IQMPO> > COp;

    for(int i = 1; i <= L; ++i) {
        std::vector<IQMPO> Ci;
        for(int j = 1; j <= L; ++j) {
            HamBuilder<IQTensor> bdb(model);
            bdb.set(model.op("Bdag",i),i,model.op("B",j),j);
            Ci.push_back(bdb);
        }
        COp.push_back(Ci);
    }

    std::cout << format("Initial energy = %.5f", psiHphi(psi,H,psi)) << std::endl;

    /*Sweeps sweeps(10);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-6,1E-7,0.0;*/
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
        n(j) = Dot(conj(primed(psi.A(j),Site)),model.op("N",j)*psi.A(j));
        n2(j) = Dot(conj(primed(psi.A(j),Site)),model.op("N2",j)*psi.A(j));
        Ntot += n(j);
        for(int k = 1; k <= L; ++k) {
            C(j, k) = psiHphi(psi,COp[j-1][k-1],psi);
        }
    }
    
    string outputfilename(argv[2]);
    ofstream outputfile(outputfilename);
    
    outputfile << "Ei ";
    for(int j = 0; j < Ei.size(); ++j) {
        outputfile << format("%.10e ", Ei[j]);
    }
    outputfile << endl;
    outputfile << "E0 " << format("%.10e", E0) << endl;
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

    cout << "Density:" << endl;
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
    }

    return 0;
}
