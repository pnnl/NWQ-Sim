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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "config.h"
#include "util.h"
#include "dmsim_nvgpu_mpi.cuh"

namespace py = pybind11;
using namespace NWQSim;
using namespace DMSim;

PYBIND11_MODULE(libdmsim, m) 
{
    py::class_<Simulation>(m, "Simulation")
        .def(py::init<>()) 
        //Basis Gate definition
        .def("X",    &Simulation::X)
        .def("ID",   &Simulation::ID)
        .def("RZ",   &Simulation::RZ)
        .def("SX",   &Simulation::SX)
        .def("CX",   &Simulation::CX)
        .def("M",    &Simulation::M)
        .def("MA",   &Simulation::MA)
        .def("RESET",&Simulation::RESET)
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
                for (unsigned i=0; i<repetition; i++) rtn.append(m_rtn[i]);
                return rtn;})
        //Get result density matrix
        .def("get_dm", [](Simulation &s) -> py::array_t<complex<ValType>, py::array::c_style | py::array::forcecast> {
                py::array_t<complex<ValType>> result = py::array_t<complex<ValType>>(s.dim);
                py::buffer_info buf = result.request();
                complex<double>* val = (complex<double>*)buf.ptr;
                for (IdxType i=0; i<s.dim; i++) 
                {
                    val[i].real(s.sv_real_cpu[i]);
                    val[i].imag(s.sv_imag_cpu[i]);
                }
                IdxType size = (IdxType)1<<(s.n_qubits);
                result.resize({size,size});
                return result;
            })
        ;
}
