#ifndef BOSE_HUBBARD_MODEL_H
#define BOSE_HUBBARD_MODEL_H

#include <model.h>

class BoseHubbardModel : public Model
{
public:

    BoseHubbardModel();

    BoseHubbardModel(int N,
                int nmax,
                const OptSet& opts = Global::opts());

    /*IQIndexVal
    Em(int i) const;

    IQIndexVal
    Occ(int i) const;

    IQIndexVal
    Dou(int i) const;

    IQIndexVal
    EmP(int i) const;

    IQIndexVal
    OccP(int i) const;

    IQIndexVal
    DouP(int i) const;*/

private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const OptSet& opts = Global::opts()) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();


    //Data members -----------------

    int N_;
    int nmax_;

    std::vector<IQIndex> site_;

};

inline BoseHubbardModel::
BoseHubbardModel()
    : N_(-1),
      nmax_(-1)
{ }

inline BoseHubbardModel::
BoseHubbardModel(int N, int nmax, const OptSet& opts)
    : N_(N),
      nmax_(nmax),
      site_(N_+1)
{
    constructSites();
}

void inline BoseHubbardModel::
constructSites()
{
    for(int j = 1; j <= N_; ++j) {
        std::vector<IndexQN> indices;
        for(int n = 0; n <= nmax_; n++) {
            indices.push_back(IndexQN(Index(nameint("n = ", n) + nameint(" for site ", j),1,Site), QN(0,n,n%2)));
        }
        site_.at(j) = IQIndex(nameint("BoseHubbard site=",j),indices);
    }
}

void inline BoseHubbardModel::
doRead(std::istream& s)
{
    s.read((char*) &N_,sizeof(N_));
    s.read((char*) &nmax_,sizeof(nmax_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j)
        site_.at(j).read(s);
}

void inline BoseHubbardModel::
doWrite(std::ostream& s) const
{
    s.write((char*) &N_,sizeof(N_));
    s.write((char*) &nmax_,sizeof(nmax_));
    for(int j = 1; j <= N_; ++j)
        site_.at(j).write(s);
}

int inline BoseHubbardModel::
getN() const
{
    return N_;
}

inline const IQIndex& BoseHubbardModel::
getSi(int i) const
{
    return site_.at(i);
}

inline IQIndexVal BoseHubbardModel::
getState(int i, const String& state) const
{
    for(int n = 0; n <= nmax_; n++) {
        if(state == nameint("",n)) {
            return getSi(i)(n+1);
        }
    }
    Error("State " + state + " not recognized");
    return IQIndexVal();
}

inline IQTensor BoseHubbardModel::
getOp(int i, const String& opname, const OptSet& opts) const
{
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    std::vector<IQIndexVal> iv, ivP;

    for(int n = 0; n <= nmax_; n++) {
        iv.push_back(s(n+1));
        ivP.push_back(sP(n+1));
    }

    IQTensor Op(conj(s),sP);

    if(opname == "N") {
        for(int n = 0; n <= nmax_; n++) {
            Op(iv[n],ivP[n]) = n;
        }
    } else    if(opname == "N2") {
        for(int n = 0; n <= nmax_; n++) {
            Op(iv[n],ivP[n]) = n*n;
        }
    } else if(opname == "OS") {
        for(int n = 0; n <= nmax_; n++) {
            Op(iv[n],ivP[n]) = n*(n-1);
        }
    }//on-site N(N-1)
    else if(opname == "B") {
        for(int n = 0; n <= nmax_-1; n++) {
            Op(iv[n+1],ivP[n]) = sqrt(n+1);
        }
    } else if(opname == "Bdag") {
        for(int n = 0; n <= nmax_-1; n++) {
            Op(iv[n],ivP[n+1]) = sqrt(n+1);
        }
    } else {
        Error("Operator " + opname + " name not recognized");
    }

    return Op;
}

#endif
