#include <numeric>
#include <cassert>
#include <string>
#include <sstream>

#include <tamm/errors.hpp>
#include <tamm/symbol.hpp>
#include <tamm/tensor.hpp>
#include <tamm/tiled_index_space.hpp>
#include <tamm/index_space.hpp>
#include <tamm/execution_context.hpp>
#include <tamm/scheduler.hpp>
#include <tamm/tamm_io.hpp>


#include <tamm/rmm_memory_pool.hpp>
#include <tamm/gpu_streams.hpp>

#include <tamm/tamm.hpp>

#include <Eigen/Dense>
#include <gsl/span>
#include <iostream>

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  GA_Initialize();

  // -----------------------------
  // 1) Initialize RMM pools
  // -----------------------------
  //  1GiB host pool, 1GiB device pool (pick sizes appropriate for your problem)

  // --------------------------------
  // Build a small 4×6 test‐tensor in TAMM
  // --------------------------------
  size_t dim1 = 4, dim2 = 6;
  tamm::IndexSpace is1({0,(tamm::Index)dim1-1}), is2({0,(tamm::Index)dim2-1});
  tamm::TiledIndexSpace tis1(is1,1), tis2(is2,1);
  tamm::Tensor<double> A({tis1,tis2}), B({tis2,tis1}), C({tis1,tis1});

  tamm::ProcGroup::self_ga_pgroup(true);
  auto pg = tamm::ProcGroup::create_world_coll();

  // Use the LOCAL memory manager, so that TAMM will use RMM+GPU
  tamm::ExecutionContext ec{pg,
                            tamm::DistributionKind::nw,
                            tamm::MemoryManagerKind::ga};

  A.allocate(&ec);
  B.allocate(&ec);
  C.allocate(&ec);

  // fill A and B on host
  size_t idx = 0;
  A.loop_nest().iterate([&](auto const& bid){
    double v = double(idx++);
    A.put(bid, gsl::span<double>(&v,1));
  });
  idx = 0;
  B.loop_nest().iterate([&](auto const& bid){
    double v = double(idx++);
    B.put(bid, gsl::span<double>(&v,1));
  });

  // --------------------------------
  // 2) Schedule a tiny GEMM on the GPU:
  //     C[i,j] = sum_k A[i,k] * B[k,j]
  // --------------------------------
  
  tamm::Scheduler sched{ec};
  sched( C("i","j") = 1.0 * A("i","k") * B("k","j"),
       "gemm",
       tamm::ExecutionHW::GPU );
  sched.execute(tamm::ExecutionHW::CPU);

  // --------------------------------
  // 3) Copy C back to Eigen on host and do your SVD
  // --------------------------------
  Eigen::MatrixXd mat(dim1,dim1);
  C.loop_nest().iterate([&](auto const& bid){
    double v;
    C.get(bid, gsl::span<double>(&v,1));
    mat(bid[0],bid[1]) = v;
  });

  // do host‐side SVD
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU|Eigen::ComputeThinV);
  std::cout << "Singular values:\n" << svd.singularValues() << "\n\n";

  // --------------------------------
  // tear down
  // --------------------------------
  C.deallocate();
  B.deallocate();
  A.deallocate();

  GA_Terminate();
  MPI_Finalize();
  return 0;
}
