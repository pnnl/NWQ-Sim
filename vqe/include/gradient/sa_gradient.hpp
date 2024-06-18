#ifndef GRADIENT
#define GRADIENT
#include <string>
#include <vector>
#include <random>


namespace NWQSim {
  using ValType = double;
  using IdxType = long long;
  namespace VQE {
  using Function = std::function<double(const std::vector<ValType>&)>; // Return a gradient (or function value) given a parameter set
  class GradientEstimator {
    protected:
      std::string name;
    public:
      GradientEstimator(std::string _name): name(_name) {};
      virtual void estimate(Function function_oracle, const std::vector<ValType>& x, std::vector<ValType>& grad, ValType epsilon) {
        throw std::runtime_error("First order oracle not implemented for this class");
      };
      virtual void estimate(Function function_oracle, const std::vector<ValType>& x, std::vector<ValType>& grad, ValType epsilon, IdxType n_evals) {
        throw std::runtime_error("First order oracle not implemented for this class");
      };
      std::string get_name() const { return name; }
  };
  class SPSA: public GradientEstimator {
    /**
     * Simultaneous Perturbation Stochastic Approximation (SPSA) Gradient Estimator
     *  Generates a random Bernoulli perturbation vector to estimate the full gradient using two function evaluations.
     * 
    */
    protected:
      std::mt19937_64 random_engine;
      std::bernoulli_distribution distribution;
    public:
      SPSA(IdxType seed): GradientEstimator("SPSA"), random_engine(seed), distribution(0.5) {}; 

      virtual void estimate(Function function_oracle, const std::vector<ValType>& x, std::vector<ValType>& grad, ValType epsilon, IdxType n_trials = 1) override {
        std::vector<ValType> delta (x.size());
        grad.resize(x.size());
        std::fill(grad.begin(), grad.end(), 0);
        for (size_t i = 0; i < n_trials; i++) {

          std::generate(delta.begin(), delta.end(), [&] () {return epsilon * (distribution(random_engine) ? 1. : -1.);});
          // Two perturbation vectors, one adding the perturbation and one subtracting
          std::vector<ValType> x1 (x.size());
          std::transform(x.begin(), x.end(), delta.begin(), x1.begin(), [&] (ValType v1, ValType v2) {return v1 + v2;} );
          std::vector<ValType> x2 (x.size());
          std::transform(x.begin(), x.end(), delta.begin(), x2.begin(), [&] (ValType v1, ValType v2) {return v1 - v2;} );

          ValType v1 = function_oracle(x1);
          ValType v2 = function_oracle(x2);
          for (size_t j = 0; j < delta.size(); j++) {
            grad[j] += (v1 - v2) /  (2 * n_trials * delta[j]);
          }
        }
      };

  };
  }; // namespace VQE
}; // namespace NWQSim




#endif