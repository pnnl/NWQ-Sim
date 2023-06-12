/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/

#include <utility>

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Optimizer.hpp"
#include "Observable.hpp"
#include "Algorithm.hpp"
#include "PauliOperator.hpp"
#include "xacc_observable.hpp"
#include "AlgorithmGradientStrategy.hpp"
#include "../nwq_accelerator.hpp"

using namespace xacc;
using namespace xacc::quantum;

#include <fstream>
#include <sstream>
#include <string>
std::string readFile(const std::string &filename)
{
  std::ifstream file(filename);
  if (!file)
  {
    throw std::runtime_error("Failed to open file");
  }

  std::ostringstream ss;
  std::string line;
  while (std::getline(file, line))
  {
    ss << line;
  }

  return ss.str();
}

int main(int argc, char **argv)
{
  xacc::Initialize(argc, argv);
  xacc::set_verbose(true);

  // Accelerator:
  // auto acc = xacc::getAccelerator("qpp");
  std::shared_ptr<xacc::Accelerator> acc = std::make_shared<xacc::quantum::NWQAccelerator>();

  auto optimizer_vqe = xacc::getOptimizer("nlopt", {std::make_pair("nlopt-optimizer", "l-bfgs")});

  auto str = readFile("be-5.txt");
  auto be_vqe = xacc::quantum::getObservable("fermion", str);

  auto tmp = xacc::getService<Instruction>("uccsd"); // std::make_shared<QFT>();
  auto uccsd = std::dynamic_pointer_cast<CompositeInstruction>(tmp);
  uccsd->expand({std::make_pair("ne", 4), std::make_pair("nq", 10),
                 std::make_pair("pool", "qubit-pool")});

  // Get the VQE Algorithm and initialize it
  auto vqe = xacc::getAlgorithm("vqe");
  vqe->initialize({std::make_pair("ansatz", uccsd),
                   std::make_pair("observable", be_vqe),
                   std::make_pair("accelerator", acc),
                   std::make_pair("optimizer", optimizer_vqe)});

  // Allocate some qubits and execute
  auto buffer = xacc::qalloc(10);
  vqe->execute(buffer);

  std::cout << buffer->getInformation("opt-val").as<double>() << std::endl;

  std::cout << "XACC VQE (BE_5 Hamiltonian) with " << acc->name() << " results:\n\n";

  auto keys = buffer->listExtraInfoKeys();

  for (auto k : keys)
  {
    std::cout << k << ": " << (*buffer)[k].toString() << "\n\n";
  }

  xacc::Finalize();
}