#ifndef BOSE_HUBBARD_HAMILTONIAN_H
#define BOSE_HUBBARD_HAMILTONIAN_H

#include <sstream>
#include <mpo.h>
#include <hambuilder.h>

#include "BoseHubbardSiteSet.h"

using namespace itensor;

class BoseHubbardHamiltonian
{
public:

    BoseHubbardHamiltonian(const BoseHubbardSiteSet& sites,
                           const OptSet& opts = Global::opts());

    std::vector<Real>
    t() const {
        return t_;
    }
    void
    t(Real val) {
        initted_ = false;
        t_.assign(sites_.N(),val);
    }
    void
    t(std::vector<Real> val) {
        initted_ = false;
        t_ = val;
    }

    std::vector<Real>
    U() const {
        return U_;
    }
    void
    U(Real val) {
        initted_ = false;
        U_.assign(sites_.N(),val);
    }
    void
    U(std::vector<Real> val) {
        initted_ = false;
        U_ = val;
    }

    std::vector<Real>
    mu() const {
        return mu_;
    }
    void
    mu(Real val) {
        initted_ = false;
        mu_.assign(sites_.N(),val);
    }
    void
    mu(std::vector<Real> val) {
        initted_ = false;
        mu_ = val;
    }

    operator MPO() {
        init_();
        return H.toMPO();
    }

    operator IQMPO() {
        init_();
        return H;
    }

private:

    ///////////////////
    //
    // Data Members

    const BoseHubbardSiteSet& sites_;
    bool initted_;
    std::vector<Real> t_,U_,mu_;
    IQMPO H;

    //
    //////////////////

    void
    init_();

}; //class BoseHubbardHamiltonian

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

double stringtodouble(std::string& s)
{
    return std::stod(s);
}

inline BoseHubbardHamiltonian::
BoseHubbardHamiltonian(const BoseHubbardSiteSet& sites,
                       const OptSet& opts)
    :
    sites_(sites),
    initted_(false),
    t_(sites.N()),
    U_(sites.N()),
    mu_(sites.N())
{
    const int Ns = sites_.N();

    std::string tstr = opts.getString("t", "1");
    std::vector<std::string> tstrs = split(tstr, ',');
    tstrs.resize(Ns, tstrs.back());
    std::transform(tstrs.begin(), tstrs.end(), t_.begin(), stringtodouble);

    std::string Ustr = opts.getString("U", "1");
    std::vector<std::string> Ustrs = split(Ustr, ',');
    Ustrs.resize(Ns, Ustrs.back());
    std::transform(Ustrs.begin(), Ustrs.end(), U_.begin(), stringtodouble);

    std::string mustr = opts.getString("mu", "0");
    std::vector<std::string> mustrs = split(mustr, ',');
    mustrs.resize(Ns, mustrs.back());
    std::transform(mustrs.begin(), mustrs.end(), mu_.begin(), stringtodouble);

}

void inline BoseHubbardHamiltonian::
init_()
{
    if(initted_) return;

    H = IQMPO(sites_);

    const int Ns = sites_.N();
    const int k = 4;
	const int nmax = sites_.nmax();

    std::vector<IQIndex> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) {
        std::vector<IndexQN> indices;
        for(int n = 0; n <= nmax; n++) {
            indices.push_back(IndexQN(Index(nameint("n = ", n) + nameint(" for site ", l),1), QN(0,n,n%2)));
        }
        links.at(l) = IQIndex(nameint("BoseHubbard site=",l),indices);
    }

    for(int n = 1; n <= Ns; ++n) {
        IQTensor& W = H.Anc(n);
        IQIndex row = dag(links[n-1]), col = links[n];

        W = IQTensor(dag(sites_.si(n)),sites_.siP(n),row,col);

        //Identity strings
        W += sites_.op("Id",n) * row(1) * col(1);
        W += sites_.op("Id",n) * row(k) * col(k);

        //Hopping terms -t*(b^d_i b_{i+1} + b_i b^d_{i+1})
        W += sites_.op("Bdag",n) * row(1) * col(2) * (-t_[n-1]);
        W += sites_.op("B",n) * row(1) * col(3) * (-t_[n-1]);
        W += sites_.op("B",n) * row(2) * col(k);
        W += sites_.op("Bdag",n) * row(3) * col(k);

        //on-site terms U/2 * n_i(n_i-1) + mu * n_i
        W += sites_.op("N",n) * row(1) * col(k) * (-mu_[n-1]);
        W += sites_.op("OS",n) * row(1) * col(k) * (0.5*U_[n-1]);//OS = N(N-1)
    }

    H.Anc(1) *= IQTensor(links.at(0)(1));
    H.Anc(Ns) *= IQTensor(dag(links.at(Ns))(k));
    
    initted_ = true;
}

#endif
