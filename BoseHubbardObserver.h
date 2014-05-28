#ifndef BOSE_HUBBARD_OBSERVER_H
#define BOSE_HUBBARD_OBSERVER_H

#include <DMRGObserver.h>

using namespace itensor;

class BoseHubbardObserver : public DMRGObserver
{
public:
    BoseHubbardObserver(const OptSet& opts = Global::opts()): DMRGObserver(opts) {}

    void virtual
    measure(int N, int sw, int ha, int b, const Spectrum& spec, Real energy,
            const OptSet& opts = Global::opts());
};

void inline BoseHubbardObserver::
measure(int N, int sw, int ha, int b, const Spectrum& spec, Real energy,
        const OptSet& opts)
{
    std::cout << "observed" << std::endl;
}

#endif
