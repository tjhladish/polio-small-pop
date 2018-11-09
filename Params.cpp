#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include "Params.h"

// struct Params {
//     double recovery;  // rate; I1->R & IR->R
//     double beta;      // infection rate per capita; S->I1, P->IR
//     double birth;     // birth rate; equivalent to death in const pop model
//     double death;     // rate; all classes -> S (via birth)
//     double kappa;     // scaling parameter for pre-exposed class
//     double rho;       // rate; waning R->P
//     std::vector<double> Population;
// };

Params parseParams(std::string jsonfile) {
  Params result = {};
  std::ifstream jsonstream(jsonfile);
  
  nlohmann::json pars;   // will contains the par value after parsing.
  jsonstream >> pars;
  
  result.recovery = pars["recovery"];
  result.beta     = pars["beta"];
  result.birth    = pars["birth"];
  result.death    = pars["death"];
  result.kappa    = pars["kappa"];
  result.rho      = pars["rho"];
  std::vector<double> tmp = pars["pop"];
  result.Population = tmp;
    // par["pops"].asFloat();
  // should be returning an Option[Params]?
  return result;
}