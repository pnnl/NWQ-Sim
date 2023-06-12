#include <utility>
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Optimizer.hpp"
#include "Observable.hpp"
#include "Algorithm.hpp"
#include "PauliOperator.hpp"
#include "xacc_observable.hpp"
#include "AlgorithmGradientStrategy.hpp"
#include "../../nwq_accelerator.hpp"

using namespace xacc;
using namespace xacc::quantum;

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

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

void processCommandLineArguments(int argc, char *argv[], std::string &ham_str)
{
  if (argc != 3)
  {
    std::cout << "Usage: program_name [-s/-d] index\n";
    exit(1);
  }

  std::string flag = argv[1];
  int index = std::stoi(argv[2]);

  // Array of filenames
  std::string filenames[] = {
      "pV5Z",
      "pV6Z",
      "pVDZ",
      "pVQZ",
      "pVTZ",
      // Add more filenames as needed
  };

  // Check the flag and retrieve the filename based on the index
  std::string filename;
  if (flag == "-s" && index >= 0 && index < sizeof(filenames) / sizeof(filenames[0]))
  {
    filename = "ham/single/" + filenames[index] + ".txt";
  }
  else if (flag == "-d" && index >= 0 && index < sizeof(filenames) / sizeof(filenames[0]))
  {
    filename = "ham/double/" + filenames[index] + ".txt";
  }
  else
  {
    std::cout << "Invalid flag or index.\n";
    exit(1);
  }

  std::ifstream file(filename);
  if (!file)
  {
    std::cout << "Failed to open file: " << filename << std::endl;
    exit(1);
  }
  else
  {
    std::cout << "Opened file: " << filename << std::endl;
    ham_str = readFile(filename);
  }
}

int main(int argc, char **argv)
{
  std::string ham_str;
  processCommandLineArguments(argc, argv, ham_str);

  xacc::Initialize(argc, argv);

  xacc::set_verbose(true);

  // Accelerator:
  std::shared_ptr<xacc::Accelerator> acc = std::make_shared<xacc::quantum::NWQAccelerator>();

  auto optimizer_vqe = xacc::getOptimizer("nlopt", {std::make_pair("nlopt-optimizer", "cobyla")});
  auto adapt_vqe = xacc::getService<xacc::Algorithm>("adapt");

  int nOrbitals = 6;
  int nElectrons = 10;
  auto pool_vqe = "qubit-pool";
  auto subAlgo_vqe = "vqe";

  auto buffer_vqe = xacc::qalloc(nOrbitals * 2);

  auto h2o_vqe = xacc::quantum::getObservable("fermion", ham_str);

  adapt_vqe->initialize({
      std::make_pair("accelerator", acc),
      std::make_pair("observable", h2o_vqe),
      std::make_pair("optimizer", optimizer_vqe),
      std::make_pair("pool", pool_vqe),
      std::make_pair("n-electrons", nElectrons),
      std::make_pair("sub-algorithm", subAlgo_vqe),
      std::make_pair("maxiter", 5),
      std::make_pair("print-operators", true),
  });

  std::cout << "Algorithm Started" << std::endl;
  adapt_vqe->execute(buffer_vqe);
  std::cout << buffer_vqe->getInformation("opt-val").as<double>() << std::endl;

  xacc::Finalize();
}