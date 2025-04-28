#include "itensor/all.h"

using namespace itensor;

int
main()
    {
    // Note: dimension comes first, then the name/tag
    Index i(2,"i"), j(2,"j");

    ITensor T(i,j);
    T.set(i=1, j=2, 3.14);

    PrintData(T);
    int N = 10;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto state = InitState(sites);
    auto B = MPS(state);

    std::cout<<B<<std::endl;

    return 0;
    }

