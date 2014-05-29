#ifndef BOSE_HUBBARD_OBSERVER_H
#define BOSE_HUBBARD_OBSERVER_H

#include <DMRGObserver.h>

using namespace itensor;

template<class Tensor>
class BoseHubbardObserver : public DMRGObserver<Tensor>
{
public:
    BoseHubbardObserver(const MPSt<Tensor>& psi, const OptSet& opts = Global::opts()): DMRGObserver<Tensor>(psi, opts), N(psi.N()) {}

    void virtual
    measure(const OptSet& opts = Global::opts());

    std::vector<Real>& getEnergies() {
        return energies_;
    }

private:
    int N;
    std::vector<Real> energies_;
};

template<class Tensor>
void inline BoseHubbardObserver<Tensor>::
measure(const OptSet& opts)
{
    const int b = opts.getInt("AtBond",1);
    const int ha = opts.getInt("HalfSweep",0);
    const Real energy = opts.getReal("Energy",0);
    if(b == N/2 && ha == 2) {
        energies_.push_back(energy);
    }
}

#endif
