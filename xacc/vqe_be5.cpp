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
#include <typeinfo>
#include <iostream>
#include <string>
#include <chrono>
#include <ctime>

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Optimizer.hpp"
#include "Observable.hpp"
#include "Algorithm.hpp"
#include "xacc_observable.hpp"
#include "AlgorithmGradientStrategy.hpp"
#include "nwq_accelerator.hpp"

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

  
  // Get reference to the Accelerator
  std::shared_ptr<xacc::Accelerator> accelerator = std::make_shared<xacc::quantum::NWQAccelerator>();
  accelerator->updateConfiguration({std::make_pair("vqe_mode", true)});

  auto str = readFile("../xacc/examples/be-5.txt");
  auto be_vqe = xacc::quantum::getObservable("fermion", str);

  auto optimizer_vqe = xacc::getOptimizer("nlopt", {std::make_pair("nlopt-optimizer", "cobyla")});

  auto adapt_vqe = xacc::getService<xacc::Algorithm>("adapt");
  int nElectrons = 4;
  int nOrbitals = 5;
  auto pool_vqe = "qubit-pool";
  auto subAlgo_vqe = "vqe";

  auto buffer_vqe = xacc::qalloc(nOrbitals * 2);

  adapt_vqe->initialize({
      std::make_pair("accelerator", accelerator),
      std::make_pair("observable", be_vqe),
      std::make_pair("optimizer", optimizer_vqe),
      std::make_pair("pool", pool_vqe),
      std::make_pair("n-electrons", nElectrons),
      std::make_pair("sub-algorithm", subAlgo_vqe),
      std::make_pair("maxiter", 20),
      std::make_pair("print-operators", true),
  });

  std::cout << "started" << std::endl;
  adapt_vqe->execute(buffer_vqe);
  std::cout << buffer_vqe->getInformation("opt-val").as<double>() << std::endl;

  xacc::Finalize();
}