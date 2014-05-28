#ifndef BOSE_HUBBARD_OBSERVER_H
#define BOSE_HUBBARD_OBSERVER_H

#include <DMRGObserver.h>

using namespace itensor;

template<class Tensor>
class BoseHubbardObserver : public DMRGObserver<Tensor>
{
public:
    BoseHubbardObserver(const MPSt<Tensor>& psi, const OptSet& opts = Global::opts()): DMRGObserver<Tensor>(psi, opts) {}

    void virtual
    measure(const OptSet& opts = Global::opts());
    
    std::vector<Real>& getEnergies() { return energies_; }
    
private:
    std::vector<Real> energies_;
};

template<class Tensor>
void inline BoseHubbardObserver<Tensor>::
measure(const OptSet& opts)
{
    const Real energy = opts.getReal("Energy",0);
    energies_.push_back(energy);
}

#endif
