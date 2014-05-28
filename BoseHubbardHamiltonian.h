#ifndef BOSE_HUBBARD_HAMILTONIAN_H
#define BOSE_HUBBARD_HAMILTONIAN_H

#include <mpo.h>
#include <hambuilder.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>


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
    std::vector<Real> t_,U_,mu_;
    bool initted_;
    IQMPO H;

    //
    //////////////////

    void
    init_();

    bool
    isReal(std::string s) {
        try         {
            boost::lexical_cast<Real>(s);
        } catch (...) {
            return false;
        }

        return true;
    }

}; //class BoseHubbardChain

inline BoseHubbardHamiltonian::
BoseHubbardHamiltonian(const BoseHubbardSiteSet& model,
                       const OptSet& opts)
    :
    model_(model),
    initted_(false)
{
    const int Ns = model_.N();

    std::string tstr = opts.getString("t","0.01");
    std::string Ustr = opts.getString("U","1");
    //std::string mustr = opts.getString("mu","0");
    std::string mustr = opts.getString("mu","0.0244067519637,0.107594683186,0.0513816880358,0.0224415914984,-0.0381726003305,0.0729470565333,-0.0312063943687,0.195886500391,0.231831380251,-0.0582792405871,0.145862519041,0.0144474598765");

    if(isReal(tstr)) {
        t_.assign(Ns, boost::lexical_cast<Real>(tstr));
    } else {
        boost::tokenizer<boost::escaped_list_separator<char> > tok(tstr);
        for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator iter = tok.begin(); iter != tok.end(); ++iter) {
            t_.push_back(boost::lexical_cast<Real>(*iter));
        }
    }
    if(isReal(Ustr)) {
        U_.assign(Ns, boost::lexical_cast<Real>(Ustr));
    } else {
        boost::tokenizer<boost::escaped_list_separator<char> > tok(Ustr);
        for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator iter = tok.begin(); iter != tok.end(); ++iter) {
            U_.push_back(boost::lexical_cast<Real>(*iter));
        }
    }
    if(isReal(mustr)) {
        mu_.assign(Ns, boost::lexical_cast<Real>(mustr));
    } else {
        boost::tokenizer<boost::escaped_list_separator<char> > tok(mustr);
        for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator iter = tok.begin(); iter != tok.end(); ++iter) {
            mu_.push_back(boost::lexical_cast<Real>(*iter));
        }
    }
}

void inline BoseHubbardHamiltonian::
init_()
{
    if(initted_) return;

    const int Ns = model_.N();

    H = IQMPO(model_);

    for(int n = 1; n <= Ns; ++n) {
        HamBuilder<IQTensor> muH(model_, model_.op("N",n), n);
        H.plusEq(-mu_[n-1] * muH);
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
