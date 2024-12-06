#include "circuit/ansatz.hpp"
#include "circuit/dynamic_ansatz.hpp"
#include "observable/fermionic_operator.hpp"
#include "utils.hpp"


#include <set>
#include <tuple>
// #include <map>




namespace NWQSim {
  namespace VQE {

    class SingletGSD: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        IdxType n_singles;
        IdxType n_doubles;
        IdxType trotter_n;
        IdxType unique_params;
        IdxType symm_level;
        Transformer qubit_transform;
        /** 
         * Enforce symmetries for each term. Each fermionic term will have one symmetry entry. If no symmetries are enforced, 
         * symmetries[i] = {{i, 1.0}}; Otherwise, symmetries[i] = {{j, 1.0}, {k, -1.0}} denotes that theta_i must be equal to theta_j - theta_k
         */
        std::vector<std::vector<std::pair<IdxType, ValType> > > symmetries;
        std::vector<IdxType> fermion_ops_to_params; // map from fermion operators to parameters (used in update)
        std::vector<std::vector<FermionOperator> > fermion_operators;
        std::set<std::tuple< IdxType, IdxType, IdxType, IdxType >> existing_tuples; // MZ: for recording symmetry
        // std::map<std::tuple<IdxType, IdxType, IdxType, IdxType>, std::vector<std::pair<IdxType, double>>> existing_terms; // MZ: for recording symmetry
        
        int counting_doubles (int n_spatial) {
          int doublt_term1_counter = 0;
          int doublt_term2_counter = 0;
          int doublt_term4_counter = 0;
          int doublt_term6_counter = 0;

          /*===========Double Excitations===========*/
          // Singlets and Triplets
          int rs = -1;
          for (IdxType r = 0; r < env.n_spatial; r++) {
            for (IdxType s = r; s < env.n_spatial; s++) {
              int pq = -1;
              for (IdxType p = 0; p < env.n_spatial; p++) {
                for (IdxType q = p; q < env.n_spatial; q++) {
                  pq += 1;
                  if (rs > pq) {
                    continue; 
                  }
                  if ( (p == r) && (q == s) ) {
                    continue;
                  }
                  if (p!=q) {
                    if (r == s) { // only singlet with two terms, not sure Group 3 or Group 4 cases
                      doublt_term2_counter += 1;
                      continue;
                    }
                    if ( ((r!=s)&&(q!=r)) || ((p==r)&&(q!=s)&(r!=s)) || ((p==s)&&(q!=r)&(r!=s)) || ((q==s)&&(p!=r)) || ((q==r)&&(p!=s))) {
                      doublt_term4_counter += 1;
                      doublt_term6_counter += 1;
                    }
                  }
                  // Group 3 to 5
                  if (p == q) {
                    // Group 3
                    if ((q != r)&&(r!=s)) {
                      // only singlet
                      doublt_term2_counter += 1;
                      continue;
                    }
                    // Group 4
                    if ((q == r)&&(r!=s)) {
                      // only singlet
                      doublt_term2_counter += 1;
                      continue;
                    }
                    // Group 5
                    if ((q != r)&&(r==s)) {
                      // only singlet
                      doublt_term1_counter += 1;
                      continue;
                    }
                  } // p == q
                } // q
              } // p
            } // s
          } // r
          return doublt_term1_counter+2*doublt_term2_counter+4*doublt_term4_counter+6*doublt_term6_counter;
          // return doublt_term1_counter+doublt_term2_counter+doublt_term4_counter+doublt_term6_counter;
        }

      public:
        SingletGSD(const MolecularEnvironment& _env, Transformer _qubit_transform, IdxType _trotter_n = 1, IdxType _symm_level = 3): 
                                  env(_env),
                                  trotter_n(_trotter_n),
                                  symm_level(_symm_level),
                                  qubit_transform(_qubit_transform),
                                  Ansatz(2 * _env.n_spatial) {
          n_singles = 2*choose2(env.n_spatial); // MZ: multiply by 2 for alpha and beta
          n_doubles = counting_doubles(env.n_spatial); // MZ: stupid way for the counting
          fermion_operators.reserve(n_singles + n_doubles);
          symmetries = std::vector<std::vector<std::pair<IdxType, ValType> > >((n_singles + n_doubles));
          fermion_ops_to_params.resize(n_doubles + n_singles);
          std::fill(fermion_ops_to_params.begin(), fermion_ops_to_params.end(), -1);
          unique_params = 0;
          existing_tuples = {};
          // existing_terms = {};
          ansatz_name = "Singlet GSD";
        };
 
        virtual std::vector<std::string> getFermionicOperatorStrings() const override {
          std::vector<std::string> result;
          result.reserve(fermion_operators.size());
          for (auto& oplist : fermion_operators) {
            std::string opstring = "";
            bool first = true;
            for (auto& op: oplist) {
              if (!first) {
                opstring = " " + opstring;
              } else {
                first = false;
              }
              opstring = op.toString(env.n_occ, env.n_virt) + opstring;
            }
            result.push_back(opstring);
          }
          return result;
        };
        

        virtual std::vector<std::pair<std::string, ValType> > getFermionicOperatorParameters() const override {
          std::vector<std::pair<std::string, ValType> > result;
          result.reserve(fermion_operators.size());
          for (size_t i = 0; i < fermion_operators.size(); i++) {
            const auto &oplist = fermion_operators.at(i);
            const std::vector<std::pair<IdxType, ValType> > &param_expr = symmetries[i];
            ValType param = 0.0;
            for (auto& i: param_expr) {
              param += i.second * theta->at(fermion_ops_to_params[i.first]);
            }
            std::string opstring = "";
            bool first = true;
            bool is_distinct = true;
            bool is_mixed = true;
            std::vector<IdxType> indices_seen(env.n_spatial, 0);
            for (auto& op: oplist) {
              if (!first) {
                opstring = " " + opstring;
              } else {
                first = false;
                is_mixed = op.getSpin() != oplist[1].getSpin();
              }
              if (indices_seen[op.getOrbitalIndex(env)]) {
                is_distinct = false;
              }
              indices_seen[op.getOrbitalIndex(env)] += 1;
              opstring = op.toString(env.n_occ, env.n_virt) + opstring;
            }
            result.push_back(std::make_pair(opstring, param));
          }
          return result;
        };
        
        const MolecularEnvironment& getEnv() const {return env;};
        virtual IdxType numParams() const override { return unique_params; };
        virtual IdxType numOps() const override { return fermion_operators.size(); };


    void add_double_excitation(FermionOperator i, FermionOperator j, FermionOperator r, FermionOperator s,  const std::vector<std::pair<IdxType, double>>& symm_expr, bool param) {

        // use this index as the unique parameter to create the symmetry
        symmetries[fermion_operators.size()] = symm_expr;
        fermion_operators.push_back({i, j, r, s});
        // record the string::parameter mapping
        if (param) {
          fermion_ops_to_params[fermion_operators.size()-1] = unique_params++;
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
        }
    }

    void add_single_excitation(FermionOperator p, FermionOperator q,  const std::vector<std::pair<IdxType, double>>& symm_expr, bool param) {

        // use this index as the unique parameter to create the symmetry
        symmetries[fermion_operators.size()] = symm_expr;
        fermion_operators.push_back({p,q});
        // record the string::parameter mapping
        if (param) {
          fermion_ops_to_params[fermion_operators.size()-1] = unique_params++;
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
        }
    }

    //---------------------------------MZ: New functions for adding excitation terms---------------------------------
    void add_single_comb(FermionOperator pa, FermionOperator qa, FermionOperator pb, FermionOperator qb, IdxType& term) {
        add_single_excitation(qa,pa,   {{term, 0.5}}, true);
        add_single_excitation(qb,pb,  {{term, 0.5}}, false);
    }

    void add_double_singlet1t(FermionOperator i, FermionOperator j, FermionOperator r, FermionOperator s, IdxType& term) {
        add_double_excitation(r,s, i,j,{{term, 1.0/sqrt(2.0)}}, true);
    }

    void add_double_singlet2tm(FermionOperator i1, FermionOperator j1, FermionOperator r1, FermionOperator s1,
                              FermionOperator i2, FermionOperator j2, FermionOperator r2, FermionOperator s2,
                              IdxType& term) {
        add_double_excitation(r1,s1,i1,j1, {{term, 0.5}}, true);
        add_double_excitation(r2,s2,i2,j2, {{term, -0.5}}, false);
    }

    void add_double_singlet2t(FermionOperator i1, FermionOperator j1, FermionOperator r1, FermionOperator s1,
                              FermionOperator i2, FermionOperator j2, FermionOperator r2, FermionOperator s2,
                              IdxType& term) {
        add_double_excitation(r1,s1,i1,j1, {{term, 0.5}}, true);
        add_double_excitation(r2,s2,i2,j2,{{term, 0.5}}, false);
    }

    void add_double_singlet4t(FermionOperator i1, FermionOperator j1, FermionOperator r1, FermionOperator s1,
                              FermionOperator i2, FermionOperator j2, FermionOperator r2, FermionOperator s2,
                              FermionOperator i3, FermionOperator j3, FermionOperator r3, FermionOperator s3,
                              FermionOperator i4, FermionOperator j4, FermionOperator r4, FermionOperator s4,
                              IdxType& term) {
        add_double_excitation(r1,s1,i1,j1, {{term, 0.5/sqrt(2.0)}}, true);
        add_double_excitation(r2,s2,i2,j2, {{term, -0.5/sqrt(2.0)}}, false);
        add_double_excitation(r3,s3,i3,j3, {{term, -0.5/sqrt(2.0)}}, false);
        add_double_excitation(r4,s4,i4,j4, {{term, 0.5/sqrt(2.0)}}, false);

    }

      void add_double_triplet(FermionOperator i1, FermionOperator j1, FermionOperator r1, FermionOperator s1,
                              FermionOperator i2, FermionOperator j2, FermionOperator r2, FermionOperator s2,
                              FermionOperator i3, FermionOperator j3, FermionOperator r3, FermionOperator s3,
                              FermionOperator i4, FermionOperator j4, FermionOperator r4, FermionOperator s4,
                              FermionOperator i5, FermionOperator j5, FermionOperator r5, FermionOperator s5,
                              FermionOperator i6, FermionOperator j6, FermionOperator r6, FermionOperator s6,
                              IdxType& term) {
        add_double_excitation(r1,s1,i1,j1, {{term, 2.0/sqrt(24.0)}}, true);
        add_double_excitation(r2,s2,i2,j2, {{term, 1.0/sqrt(24.0)}}, false);
        add_double_excitation(r3,s3,i3,j3, {{term, 1.0/sqrt(24.0)}}, false);
        add_double_excitation(r4,s4,i4,j4, {{term, 1.0/sqrt(24.0)}}, false);
        add_double_excitation(r5,s5,i5,j5, {{term, 1.0/sqrt(24.0)}}, false);
        add_double_excitation(r6,s6,i6,j6, {{term, 2.0/sqrt(24.0)}}, false);
    }
    //------------------------------------------------------------------
    // MZ: Show this orbital should be Occupied or Virtual based on the spatial index
    auto OccVir(IdxType sp, IdxType n_occ) {
      if (sp < n_occ) { // assume first n_occ # of orbitals are occupied
        return Occupied;
      }
      return Virtual; 
    }

    // MZ: From spatial index to the index used for FermionicOperator
    IdxType InCo(IdxType sp, IdxType n_occ) {
      if (sp < n_occ) { // assume first n_occ # of orbitals are occupied
        return sp;
      }
      return sp - n_occ;
    }

   /**
    * @brief  Generate Fermionic operators for UCCSD
    * @note   Symmetry-linked operators (e.g. by anticommutation, spin reversal) share parameters
    * @retval None
    */
    void  getFermionOps() {
      int single_counter = 0;
      int double_counter = 0;
      int doublt_term1_counter = 0;
      int doublt_term2_counter = 0;
      int doublt_term4_counter = 0;
      int doublt_term6_counter = 0;
      /*===========Single Excitations===========*/
      for (IdxType p = 0; p < env.n_spatial; p++) {
        IdxType pi = InCo(p, env.n_occ);
        auto pov = OccVir(p, env.n_occ);
        FermionOperator virtual_creation_up (pi, pov, Up, Creation, env.xacc_scheme);
        FermionOperator virtual_creation_down (pi, pov, Down, Creation, env.xacc_scheme);
        for (IdxType q = p+1; q < env.n_spatial; q++) {
          IdxType qi = InCo(q, env.n_occ);
          auto qov = OccVir(q, env.n_occ);
          // creation operator
          FermionOperator occupied_annihilation_up (qi, qov, Up, Annihilation, env.xacc_scheme);
          FermionOperator occupied_annihilation_down (qi, qov, Down, Annihilation, env.xacc_scheme);
          // Alpha+beta
          IdxType term_single = fermion_operators.size();
          add_single_comb(virtual_creation_up, occupied_annihilation_up, virtual_creation_down, occupied_annihilation_down, term_single);
          single_counter += 1;
        }
      }
      /*===========Double Excitations===========*/
      // Singlets and Triplets
      int rs = -1;
      for (IdxType r = 0; r < env.n_spatial; r++) {
        IdxType ri = InCo(r, env.n_occ);
        auto rov = OccVir(r, env.n_occ);
        FermionOperator ra (ri, rov, Up, Annihilation, env.xacc_scheme);
        FermionOperator rb (ri, rov, Down, Annihilation, env.xacc_scheme);
        for (IdxType s = r; s < env.n_spatial; s++) {
          IdxType si = InCo(s, env.n_occ);
          auto sov = OccVir(s, env.n_occ);
          FermionOperator sa (si, sov, Up, Annihilation, env.xacc_scheme);
          FermionOperator sb (si, sov, Down, Annihilation, env.xacc_scheme);
          rs += 1;
          int pq = -1;
          for (IdxType p = 0; p < env.n_spatial; p++) {
            IdxType pi = InCo(p, env.n_occ);
            auto pov = OccVir(p, env.n_occ);
            FermionOperator pa (pi, pov, Up, Creation, env.xacc_scheme);
            FermionOperator pb (pi, pov, Down, Creation, env.xacc_scheme);
            for (IdxType q = p; q < env.n_spatial; q++) {
              IdxType qi = InCo(q, env.n_occ);
              auto qov = OccVir(q, env.n_occ);
              FermionOperator qa (qi, qov, Up, Creation, env.xacc_scheme);
              FermionOperator qb (qi, qov, Down, Creation, env.xacc_scheme);
              pq += 1;
              if (rs > pq) {
                continue; 
              }
              if ( (p == r) && (q == s) ) {
                continue;
              }
              if (p!=q) {
                if (r == s) { // only singlet with two terms, not sure Group 3 or Group 4 cases
                  IdxType term1 = fermion_operators.size(); // May need +1
                  add_double_singlet2tm(qb, pa, rb, ra, qa, pb, rb, ra,  term1);
                  double_counter += 1;
                  doublt_term2_counter += 1;
                  continue;
                }
                if ( ((r!=s)&&(q!=r)) || ((p==r)&&(q!=s)&(r!=s)) || ((p==s)&&(q!=r)&(r!=s)) || ((q==s)&&(p!=r)) || ((q==r)&&(p!=s))) {
                  // Singlet
                  IdxType term1 = fermion_operators.size();
                  add_double_triplet(qa,pa,sa,ra,
                                     qb,pa,sb,ra,
                                     qb,pa,sa,rb,
                                     qa,pb,sb,ra,
                                     qa,pb,sa,rb,
                                     qb,pb,sb,rb,
                                     term1);
                  double_counter += 1;
                  doublt_term6_counter += 1;
                  continue;
                }
              }
              // Group 3 to 5
              if (p == q) {
                // Group 3
                if ((q != r)&&(r!=s)) {
                  // only singlet
                  IdxType term1 = fermion_operators.size();
                  add_double_singlet2t(pb,pa,sb,ra,   pb,pa,sa,rb, term1);
                  double_counter += 1;
                  doublt_term2_counter += 1;
                  continue;
                }
                // Group 4
                if ((q == r)&&(r!=s)) {
                  // only singlet
                  IdxType term1 = fermion_operators.size();
                  add_double_singlet2t(pa,pb,sb,ra,   pb,pa,sa,rb, term1);
                  double_counter += 1;
                  doublt_term2_counter += 1;
                  continue;
                }
                // Group 5
                if ((q != r)&&(r==s)) {
                  // only singlet
                  IdxType term1 = fermion_operators.size();
                  add_double_singlet1t(pb,pa,rb,ra, term1);
                  double_counter += 1;
                  doublt_term1_counter += 1;
                  continue;
                }
              } // p == q

            } // q
          } // p
        } // s
      } // r

      int rs_t = -1;
      for (IdxType r = 0; r < env.n_spatial; r++) {
        IdxType ri = InCo(r, env.n_occ);
        auto rov = OccVir(r, env.n_occ);
        FermionOperator ra (ri, rov, Up, Annihilation, env.xacc_scheme);
        FermionOperator rb (ri, rov, Down, Annihilation, env.xacc_scheme);
        for (IdxType s = r; s < env.n_spatial; s++) {
          IdxType si = InCo(s, env.n_occ);
          auto sov = OccVir(s, env.n_occ);
          FermionOperator sa (si, sov, Up, Annihilation, env.xacc_scheme);
          FermionOperator sb (si, sov, Down, Annihilation, env.xacc_scheme);
          rs_t += 1;
          int pq_t = -1;
          for (IdxType p = 0; p < env.n_spatial; p++) {
            IdxType pi = InCo(p, env.n_occ);
            auto pov = OccVir(p, env.n_occ);
            FermionOperator pa (pi, pov, Up, Creation, env.xacc_scheme);
            FermionOperator pb (pi, pov, Down, Creation, env.xacc_scheme);
            for (IdxType q = p; q < env.n_spatial; q++) {
              IdxType qi = InCo(q, env.n_occ);
              auto qov = OccVir(q, env.n_occ);
              FermionOperator qa (qi, qov, Up, Creation, env.xacc_scheme);
              FermionOperator qb (qi, qov, Down, Creation, env.xacc_scheme);
              pq_t += 1;
              if (rs_t > pq_t) {
                continue; 
              }
              if ( (p == r) && (q == s) ) {
                continue;
              }
              if ((p == q) || (r==s)) {
                continue;
              } // singlet cases

              if ( ((r!=s)&&(q!=r)) || ((p==r)&&(q!=s)&(r!=s)) || ((p==s)&&(q!=r)&(r!=s)) || ((q==s)&&(p!=r)) || ((q==r)&&(p!=s))) {
                IdxType term2 = fermion_operators.size();
                add_double_singlet4t(qb,pa,sb,ra,
                                      qb,pa,sa,rb,
                                      qa,pb,sb,ra,
                                      qa,pb,sa,rb,
                                      term2);
                double_counter += 1;
                doublt_term4_counter += 1;
                // std::cout << "Triplet&Singlets: " << p << q << r << s << "\n" << std::endl;
                continue;
              }

            } // q
          } // p
        } // s
      } // r


      #ifndef NDEBUG
        std::cout << "Operator Stats" << std::endl;
        std::cout << "Single Excitations: " << single_counter << std::endl;
        std::cout << "Double Excitations: " << double_counter << std::endl;
        std::cout << "  -Term 1: " << doublt_term1_counter << std::endl;
        std::cout << "  -Term 2: " << doublt_term2_counter << std::endl;
        std::cout << "  -Term 4: " << doublt_term4_counter << std::endl;
        std::cout << "  -Term 6: " << doublt_term6_counter << std::endl;
      #endif

      #ifndef NDEBUG
        // Print out the generated Fermionic operators
        for (size_t i = 0; i < fermion_operators.size(); i++) {
          std::string fermi_op = "";
          for (auto i: fermion_operators[i]) {
            fermi_op += i.toString(env.n_occ, env.n_virt) + " ";
          }
          std::cout << i << ": " << fermi_op << " || [";
          for (auto sym: symmetries[i]) {
            std::cout << "{" <<sym.first << ", " << fermion_ops_to_params[sym.first] << ", " << sym.second << "}" << ", ";
          }
          std::cout << "]  " << fermion_ops_to_params[i] << std::endl;
          // if (symmetries[i].size())
          // std::cout  << std::endl;
        }
        std::cout << fermion_operators.size() << " " << unique_params << std::endl;
      #endif

    } // getFermionOps


    void  buildAnsatz()  override {
      getFermionOps();
      // assert((n_doubles + n_singles) == fermion_operators.size());
      std::cout << "Number of single excitations: " << n_singles << ", number of double excitations: " << n_doubles << std::endl;
      std::cout << "Generated " << fermion_operators.size() << " operators." << std::endl;
      theta->resize(unique_params * trotter_n);
      // exit(0);
      std::vector<std::vector<PauliOperator> > pauli_oplist;
      if (env.xacc_scheme) {
        for (IdxType i = 0; i < env.n_occ; i++) {
          X(i);
          X(i+env.n_spatial);
        }
      } else {
        for (IdxType i = 0; i < 2 * env.n_occ; i++) {
          X(i);
        }
      }
      pauli_oplist.reserve(4 * n_singles + 16 * n_doubles);
      qubit_transform(env, fermion_operators, pauli_oplist, true);  
      IdxType index = 0; // parameter index, shares parameters for Pauli evolution gates corresponding to the same Fermionic operator within the same Trotter step
      for (auto& fermionic_group: pauli_oplist) {
        bool wasused = 0;
        for (auto& pauli: fermionic_group) {
          assert (pauli.getCoeff().imag() == 0.0);
          double coeff = pauli.getCoeff().real(); 
          
          std::vector<std::pair<IdxType, ValType> > idxvals(symmetries[index].size());

          std::transform(symmetries[index].begin(), symmetries[index].end(),idxvals.begin(), 
            [&] (std::pair<IdxType, ValType> val) {
              return std::make_pair(fermion_ops_to_params[val.first], val.second);
            } );

          if (pauli.isNonTrivial() && abs(coeff) > 1e-10) { 
            ExponentialGate(pauli, OP::RZ, idxvals, 2 * coeff);
            wasused = 1;
          }
        }
        if(!wasused) {
          printf("%lld parameter not used for operator\n", index);
        }
        index++;
      }
      for (IdxType i = 0; i < trotter_n - 1; i++) {
        for (auto& fermionic_group: pauli_oplist) {
          for (auto& pauli: fermionic_group) {
            double coeff = pauli.getCoeff().real();
            if (pauli.isNonTrivial() && abs(coeff) > 1e-10)  {  
            std::vector<std::pair<IdxType, ValType> > idxvals = symmetries[index];
            std::transform(symmetries[index].begin(), symmetries[index].end(),idxvals.begin(), 
            [&, i] (std::pair<IdxType, ValType> val) {
              return std::make_pair(fermion_ops_to_params[val.first] + (i + 1) * unique_params, val.second);
            } );
            ExponentialGate(pauli, OP::RZ, idxvals, 2 * coeff);
            }
          }
          index++;
        }
      }
    }


    }; // class Singlet GSD
  
  
  };// namespace vqe
};// namespace nwqsim