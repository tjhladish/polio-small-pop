#ifndef STATETYPE_H
#define STATETYPE_H

#include <vector>
#include <string>
#include "Params.h"
#include "States.h"

enum EventType {
  FIRST_INFECTION,
  REINFECTION,
  RECOVER_FIRST_INFECTION,
  RECOVER_REINFECTION,
  WANING,
  BIRTHDEATH,
  NUM_OF_EVENT_TYPES
};

const std::vector<const EventType> eventref = {
  FIRST_INFECTION,
  REINFECTION,
  RECOVER_FIRST_INFECTION,
  RECOVER_REINFECTION,
  WANING,
  BIRTHDEATH
};

const std::vector<const std::string> eventstr = {
  "1INF",
  "REINF",
  "1REC",
  "REREC",
  "WANING",
  "MORT"
};

std::vector<double> set_rates(std::vector<int> &states, Params * p);
void update(std::vector<double> &states, Params p, EventType ev, std::vector<double> &rates, gsl_rng *rng);

#endif