// minimal_tamm_alloc_bug.cpp
#include <tamm/tamm.hpp>
#include <complex>
#include <iostream>

int main(int argc, char** argv) {
    // 1) initialize Tamm (+ MPI, GA)
    tamm::initialize(argc, argv);

    // 2) single procgroup of size=world
    auto pg = tamm::ProcGroup::create_world_coll();

    // 3) dense distribution + GA memory manager
    tamm::ExecutionContext ec{pg,
                              tamm::DistributionKind::dense,
                              tamm::MemoryManagerKind::ga};

    // 4) hard‐code three tiled index spaces of sizes 1,2,10
    tamm::IndexSpace is1{tamm::range(0,1)};   // 1 index: 0
    tamm::IndexSpace is2{tamm::range(0,2)};   // 2 indices: 0,1
    tamm::IndexSpace is3{tamm::range(0,10)};  // 10 indices: 0..9

    tamm::TiledIndexSpace tis1{is1, 1};
    tamm::TiledIndexSpace tis2{is2, 1};
    tamm::TiledIndexSpace tis3{is3, 1};

    // 5) create (1×2×10) tensor of complex<double>
    using Cplx = std::complex<double>;
    tamm::Tensor<Cplx> T{tis1, tis2, tis3};

    std::cerr << "About to allocate a (1,2,10) tensor...\n";

    // <-- this line is where you will reproduce the EXPECTS(...) failure
    T.set_dense();
    T.allocate(&ec);

    std::cerr << "Allocation succeeded (unexpected!)\n";

    // 6) tear down
    tamm::finalize();
    return 0;
}

