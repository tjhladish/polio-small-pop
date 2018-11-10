#ifndef STATETYPE_H
#define STATETYPE_H

#include <vector>
#include <string>
#include "Params.h"

enum StateType {
  S_STATE,
  I1_STATE,
  R_STATE,
  P_STATE,
  IR_STATE,
  NUM_OF_STATE_TYPES
};

const std::vector<const StateType> stateref = { S_STATE, I1_STATE, R_STATE, P_STATE, IR_STATE };

const std::vector<const std::string> statestr = {"S", "I1", "R", "P", "IR"};

std::vector<double> equilibrium_fraction(Params p);

void printResults(std::vector<double> res);

#endif