// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/dm-sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: dmsim_python.cpp
// Python wrapper via Pybind11 for DMSim
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <string>
#include <bitset>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
//#include <mpi4py/mpi4py.h>

#include "config.h"
#include "util.h"
#include "svsim_nvgpu_sin.cuh"

namespace py = pybind11;
using namespace NWQSim;

//Need to convert from py::object to MPI_Comm, and use it to init the Simulation object
//use method here: https://stackoverflow.com/questions/62309487/pybind11-init-with-lambda

// Return a MPI communicator from mpi4py communicator object
//MPI_Comm *get_mpi_comm(py::object py_comm) 
//{
//auto comm_ptr = PyMPIComm_Get(py_comm.ptr());
//if (!comm_ptr) throw py::error_already_set();
//return comm_ptr;
//}


//class Comm
//{
//private: 
//Comm(int); //private constructor
//public:
////Factory function - returned by value:
//static Comm create(int a) { return Comm(a);}
//Comm(int, int);
//}


//.def(py::init(  [](int a, int b){return Comm(a,b);}      )) 


PYBIND11_MODULE(libsvsim, m) 
{
    // initialize mpi4py's C-API
    //if (import_mpi4py() < 0) throw py::error_already_set();

    //MPI_Init(NULL, NULL);


     py::class_<Simulation>(m, "Simulation")
        .def(py::init<IdxType>()) 
        //Basis Gate definition
        .def("X",    &Simulation::X)
        .def("Y",    &Simulation::Y)
        .def("Z",    &Simulation::Z)
        .def("H",    &Simulation::H)
        .def("S",    &Simulation::S)
        .def("SDG",  &Simulation::SDG)
        .def("T",    &Simulation::T)
        .def("TDG",  &Simulation::TDG)
        .def("RI",   &Simulation::RI)
        .def("RX",   &Simulation::RX)
        .def("RY",   &Simulation::RY)
        .def("RZ",   &Simulation::RZ)
        .def("SX",   &Simulation::SX)
        .def("P",    &Simulation::P)
        .def("U",    &Simulation::U)
        .def("CX",   &Simulation::CX)
        .def("CY",   &Simulation::CY)
        .def("CZ",   &Simulation::CZ)
        .def("CH",   &Simulation::CH)
        .def("CS",   &Simulation::CS)
        .def("CSDG", &Simulation::CSDG)
        .def("CT",   &Simulation::CT)
        .def("CTDG", &Simulation::CTDG)
        .def("CRX",  &Simulation::CRX)
        .def("CRY",  &Simulation::CRY)
        .def("CRZ",  &Simulation::CRZ)
        .def("CSX",  &Simulation::CSX)
        .def("CP",   &Simulation::CP)
        .def("CU",   &Simulation::CU)
        .def("ID",   &Simulation::ID)
        .def("SWAP", &Simulation::SWAP)
        .def("M",    &Simulation::M)
        .def("MA",   &Simulation::MA)
        .def("RESET",&Simulation::RESET)
        //Composition Gates
        .def("U1",   &Simulation::U1)
        .def("U2",   &Simulation::U2)
        .def("U3",   &Simulation::U3)
        .def("CCX",  &Simulation::CCX)
        .def("CSWAP",&Simulation::CSWAP)
        .def("RXX",  &Simulation::RXX)
        .def("RYY",  &Simulation::RYY)
        .def("RZZ",  &Simulation::RZZ)

        //Simulator operations
        .def("set_seed", &Simulation::set_seed)
        .def("reset_sim", &Simulation::reset_sim)
        .def("reset_circuit", &Simulation::reset_circuit)
        .def("get_qubits", &Simulation::get_n_qubits)
        .def("get_gates", &Simulation::get_n_gates)
        .def("run", &Simulation::sim)
        .def("measure", &Simulation::measure)
        .def("measure_all",[](Simulation &s, unsigned repetition) -> py::list{
                IdxType* m_rtn = s.measure_all(repetition);
                py::list rtn;
                for (unsigned i=0; i<repetition; i++)
                {
                std::string s = std::bitset<32>(m_rtn[i]).to_string();
                rtn.append(s);
                //rtn.append(m_rtn[i]);
                }
                return rtn;})
        //This is to set comm when using python as shell for MPI calling
        //.def("set_comm", [](Simulation &s, py::object py_comm) {
        //MPI_Comm* comm = get_mpi_comm(py_comm);
        //s.set_comm(*comm);
        //})
        //Get result state vector in numpy format
        .def("get_sv", [](Simulation &s) -> py::array_t<complex<ValType>, py::array::c_style | py::array::forcecast> {
                py::array_t<complex<ValType>> result = py::array_t<complex<ValType>>(s.dim);
                py::buffer_info buf = result.request();
                complex<double>* val = (complex<double>*)buf.ptr;
                for (IdxType i=0; i<s.dim; i++) 
                {
                    val[i].real(s.sv_real_cpu[i]);
                    val[i].imag(s.sv_imag_cpu[i]);
                }
                IdxType size = (IdxType)1<<(s.n_qubits);
                result.resize({size});
                return result;
            })
        ;
}
