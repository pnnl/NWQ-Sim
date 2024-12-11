# Hamiltonians

* `H415.hamil`: linear H4 (4 spartial orbitals, 4 particles), R(H-H) = 1.5 Angstrom
* `LiH34.hamil`: LiH (6 spartial orbitals, 4 particles), R(Li-H) = 3.4 Angstrom
* `H618.hamil`: linear H6 (6 spartial orbitals, 6 particles), R(H-H) = 1.8 Angstrom
* `H622.hamil`: linear H6 (6 spartial orbitals, 6 particles), R(H-H) = 2.2 Angstrom

|  Source/Method   |      Ansatz  |     H4 (1.5 Å)    |    LiH (3.4 Å)    |     H6 (1.8 Å)    |     H6 (2.2 Å)    |
|:----------------:|:------------:|:-----------------:|:-----------------:|:-----------------:|:-----------------:|
|Eigendecomposition|  FCI energy  | -1.99615032551880 | -7.78949878539877 | -2.88756937643578 | -2.82434165191321 |
|          VQE     | UCCSD Symm 2 | -1.99476397266065 | -7.78937734881151 |                   | -2.80862380422462 |
|          VQE     |  Singlet GSD | -1.99614746813286 | -7.78948638738587 |                   | -2.82346354815650 |
|    ADAPT-VQE     |     UCCSD    | -1.99471374600654 | -7.78937885628381 |                   | -2.82356033439700 |
|    ADAPT-VQE     |  Singlet GSD | -1.99614808581472 | -7.78948638738587 |                   |                   |
| [Data from paper](10.1021/acs.jctc.9b01083) | FCI energy    |         /         | -7.78949878516939 | -2.88756937438597 | -2.82434165111152 |
| [Data from paper](10.1021/acs.jctc.9b01083)  |     UCCSD    |         /         | -7.78938765433587 | -2.87774059976972 | -2.81118282792627 |


## Example commands used for the table

VQE with UCCSD and Symmetry level 2 (0 lowest, 2 highest)
```shell
mpirun -n 4 ./vqe/nwq_vqe -f ./LiH34.hamil -p 4 -b MPI -o LN_NEWUOA -v --abstol 1e-10 --maxeval 4000 --sym 2
```

VQE with Singlet and Triplet excitations
```shell
mpirun -n 4 ./vqe/nwq_vqe -f ./LiH34.hamil -p 4 -b MPI -o LN_NEWUOA -v --abstol 1e-10 -lb -0.9 -ub 0.9 --maxeval 20000 --gsd
```

ADAPT-VQE with UCCSD ansatz pool
```shell
mpirun -n 4 ./vqe/nwq_vqe -f ./LiH34.hamil -p 4 -b MPI --abstol 1e-10 --maxeval 4000 -lb -0.9 -ub 0.9 --adapt --adapt-fvaltol -1
```

ADAPT-VQE with singlet and triplet ansatz pool (this pool currently is not normalized, unlike its VQE/QFlow version, seems like non-normalized one gives $10\%$ smaller error only for ADAPT-VQE in this example)
```shell
mpirun -n 4 ./vqe/nwq_vqe -f ./LiH34.hamil -p 4 -b MPI --abstol 1e-10 --maxeval 16000 --gsd --adapt --adapt-fvaltol -1
```


