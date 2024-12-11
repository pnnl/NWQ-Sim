# ADAPT-VQE

Note that `state->optimize()` (the VQE subroutine) is override in `vqe/include/vqe_state.hpp`, `vqe/include/svsim_vqe/sv_mpi_vqe.hpp`, and `vqe/include/svsim_vqe/sv_cuda_mpi_vqe.cuh`

Current termination flags for ADAPT-VQE:            
```shell
0 -> Reach the gradient norm tolerance
1 -> Reach the function tolerance
2 -> Reach the maximum number of ADAPT iterations
9 -> ADAPT iteration is not run successfully
```