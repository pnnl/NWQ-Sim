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
// File: python_interface.cpp
// Python wrapper via Pybind11 for NWQSim
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "config.h"
#include "util.h"

#ifdef USE_NVGPU //NVGPU
#ifdef USE_OMP//OMP
#include "dmsim_nvgpu_omp.cuh"
#elif defined USE_MPI//MPI
#include "dmsim_nvgpu_mpi.cuh"
#else //Single
#include "dmsim_nvgpu_sin.cuh"
#endif
#elif defined USE_AMDGPU //AMDGPU
#ifdef USE_OMP//OMP
#include "dmsim_amdgpu_omp.hpp"
#elif defined USE_MPI//MPI
#include "dmsim_amdgpu_mpi.hpp"
#else //Single
#include "dmsim_amdgpu_sin.hpp"
#endif
#else //CPU
#ifdef USE_OMP//OMP
#include "dmsim_cpu_omp.hpp"
#elif defined USE_MPI//MPI
#include "dmsim_cpu_mpi.hpp"
#else //Single
#include "dmsim_cpu_sin.hpp"
#endif
#endif

namespace py = pybind11;
using namespace DMSim;

PYBIND11_MODULE(libdmsim, m) 
{
    py::class_<Simulation>(m, "Simulation")
        .def(py::init<IdxType,IdxType>()) //n_qubits, n_cpus
        //Gate definition
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
        .def("CRI",  &Simulation::CRI)
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
        .def("MZ",   &Simulation::MZ)
        .def("RESET",&Simulation::RESET)
        .def("CCX",  &Simulation::CCX)
        .def("CSWAP",&Simulation::CSWAP)
        .def("U1",   &Simulation::U1)
        .def("U2",   &Simulation::U2)
        .def("U3",   &Simulation::U3)
        //Simulator operations
        .def("set_seed", &Simulation::set_seed)
        .def("reset_sim", &Simulation::reset_sim)
        .def("reset_circuit", &Simulation::reset_circuit)
        .def("get_qubits", &Simulation::get_n_qubits)
        .def("get_gates", &Simulation::get_n_gates)
        .def("run", &Simulation::sim)
        .def("measure", &Simulation::measure)
        .def("measureZ", &Simulation::measureZ)
        .def("measure_all",[](Simulation &s, unsigned repetition) -> py::list{
                IdxType* m_rtn = s.measure_all(repetition);
                py::list rtn;
                for (unsigned i=0; i<repetition; i++) rtn.append(m_rtn[i]);
                return rtn;})
        //Set noise superoperator
        .def("set_noise_1q_sop", [](Simulation &s, py::array_t<complex<ValType>, py::array::c_style | py::array::forcecast> array) -> void{
                SAFE_ALOC_HOST(s.noise_1q_real, 16*sizeof(ValType));
                SAFE_ALOC_HOST(s.noise_1q_imag, 16*sizeof(ValType));
                py::buffer_info buf = array.request();
                complex<double>* val = (complex<double>*)buf.ptr;
                for (unsigned i=0; i<16; i++) 
                {
                    s.noise_1q_real[i] = val[i].real();
                    s.noise_1q_imag[i] = val[i].imag();
                    //printf("real:%lf,imag:%lf\n",val[i].real(),val[i].imag());
                }
            })
        .def("set_noise_2q_sop", [](Simulation &s, py::array_t<complex<ValType>, py::array::c_style | py::array::forcecast> array) -> void{
                SAFE_ALOC_HOST(s.noise_2q_real, 256*sizeof(ValType));
                SAFE_ALOC_HOST(s.noise_2q_imag, 256*sizeof(ValType));
                py::buffer_info buf = array.request();
                complex<double>* val = (complex<double>*)buf.ptr;
                for (unsigned i=0; i<256; i++) 
                {
                    s.noise_2q_real[i] = val[i].real();
                    s.noise_2q_imag[i] = val[i].imag();
                }
            })
        //Get result density matrix
        .def("get_dm", [](Simulation &s) -> py::array_t<complex<ValType>, py::array::c_style | py::array::forcecast> {
                py::array_t<complex<ValType>> result = py::array_t<complex<ValType>>(s.dim);
                py::buffer_info buf = result.request();
                complex<double>* val = (complex<double>*)buf.ptr;
                for (IdxType i=0; i<s.dim; i++) 
                {
                    val[i].real(s.sv_real[i]);
                    val[i].imag(s.sv_imag[i]);
                }
                IdxType size = (IdxType)1<<(s.n_qubits);
                result.resize({size,size});
                return result;
            })
        ;
}
