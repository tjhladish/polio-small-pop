#ifndef STATETYPE_H
#define STATETYPE_H

#include <vector>
#include <string>

enum StateType {
  S_STATE,
  I1_STATE,
  R_STATE,
  P_STATE,
  IR_STATE,
  NUM_OF_STATE_TYPES
};

const std::vector<std::string> statestr = {"S", "I1", "R", "P", "IR"};

#endif