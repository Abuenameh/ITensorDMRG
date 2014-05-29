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

    BoseHubbardHamiltonian(const BoseHubbardSiteSet& model,
                           const OptSet& opts = Global::opts());

    std::vector<Real>
    t() const {
        return t_;
    }
    void
    t(Real val) {
        initted_ = false;
        t_.assign(model_.N(),val);
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
        U_.assign(model_.N(),val);
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
        mu_.assign(model_.N(),val);
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
        return H;//.toIQMPO();
    }

private:

    ///////////////////
    //
    // Data Members

    const BoseHubbardSiteSet& model_;
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
BoseHubbardHamiltonian(const BoseHubbardSiteSet& model,
                       const OptSet& opts)
    :
    model_(model),
    initted_(false),
    t_(model.N()),
    U_(model.N()),
    mu_(model.N())
{
    const int Ns = model_.N();

    std::string tstr = opts.getString("t");
    std::vector<std::string> tstrs = split(tstr, ',');
    tstrs.resize(Ns, tstrs.back());
    std::transform(tstrs.begin(), tstrs.end(), t_.begin(), stringtodouble);

    std::string Ustr = opts.getString("U");
    std::vector<std::string> Ustrs = split(Ustr, ',');
    Ustrs.resize(Ns, Ustrs.back());
    std::transform(Ustrs.begin(), Ustrs.end(), U_.begin(), stringtodouble);

    std::string mustr = opts.getString("mu");
    std::vector<std::string> mustrs = split(mustr, ',');
    mustrs.resize(Ns, mustrs.back());
    std::transform(mustrs.begin(), mustrs.end(), mu_.begin(), stringtodouble);

}

void inline BoseHubbardHamiltonian::
init_()
{
    if(initted_) return;

    const int Ns = model_.N();

    H = IQMPO(model_);

    for(int n = 1; n <= Ns; ++n) {
        HamBuilder<IQTensor> muH(model_, model_.op("N",n), n);
        if(mu_[n-1] != 0) {
            H.plusEq(-mu_[n-1] * muH);
        }
        HamBuilder<IQTensor> UH(model_, model_.op("OS",n), n);
        H.plusEq(0.5*U_[n-1] * UH);
        if(n < Ns) {
            HamBuilder<IQTensor> tH1(model_, model_.op("Bdag",n), n, model_.op("B",n+1), n+1);
            HamBuilder<IQTensor> tH2(model_, model_.op("B",n), n, model_.op("Bdag",n+1), n+1);
            H.plusEq(-t_[n-1] * tH1);
            H.plusEq(-t_[n-1] * tH2);
        }
    }

    initted_ = true;
}

#endif
