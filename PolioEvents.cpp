#include <vector>
#include <string>
#include "Params.h"
#include "States.h"

// enum EventType {
//   FIRST_INFECTION,
//   REINFECTION,
//   RECOVER_FIRST_INFECTION,
//   RECOVER_REINFECTION,
//   WANING,
//   BIRTHDEATH,
//   NUM_OF_EVENT_TYPES
// };

inline StateType from(EventType ev) {
  assert(ev != NUM_OF_EVENT_TYPES);
  assert(ev != BIRTHDEATH);
  switch (ev) {
    case FIRST_INFECTION:         return S_STATE;
    case REINFECTION:             return P_STATE;
    case RECOVER_FIRST_INFECTION: return I1_STATE;
    case RECOVER_REINFECTION:     return IR_STATE;
    case WANING:                  return R_STATE;
    default:                      exit(0);
  }
}

inline StateType to(EventType ev) {
  assert(ev != NUM_OF_EVENT_TYPES);
  switch (ev) {
    case FIRST_INFECTION:         return I1_STATE;
    case REINFECTION:             return IR_STATE;
    case RECOVER_FIRST_INFECTION: return R_STATE;
    case RECOVER_REINFECTION:     return R_STATE;
    case WANING:                  return P_STATE;
    case BIRTHDEATH:              return S_STATE;
    default: exit(0);
  }
}

inline double rate(EventType ev, std::vector<double> &states, Params * p, double foi) {
  switch (ev) {
    case FIRST_INFECTION:         return foi*states[S_STATE];
    case REINFECTION:             return foi*p->kappa*states[P_STATE];
    case RECOVER_FIRST_INFECTION: return p->recovery*states[I1_STATE];
    case RECOVER_REINFECTION:     return (p->recovery/p->kappa)*states[IR_STATE];
    case WANING:                  return p->rho*states[R_STATE];
    case BIRTHDEATH:              return p->death*(TOT - states[S_STATE]);
  }
}

void update_rates(std::vector<int> &states, Params * p, std::vector<double> &rates, std::vector<EventType> &events) {  
  double foi = p->beta/TOT*(states[I1_STATE] + p->kappa*states[IR_STATE]);
  for (auto ev : events) rates[ev] = rate(ev, states, p, foi);
};

std::vector<double> set_rates(std::vector<int> &states, Params * p) {  
  vector<double> res(NUM_OF_EVENT_TYPES, 0.0);
  update_rates(states, p, res, eventref);
  return res;
};

inline void fromto(StateType from, StateType to, std::vector<double> &states) {
  states[from]--;  states[to]++;
}

std::vector<EventType> birthdeath_rate_changes(StateType from) {
  // updates = { BIRTHDEATH, FIRST_INFECTION }; // TODO change below to pushbacks of the differing updates
  switch(from) {
    case I1_STATE:
      return { BIRTHDEATH, FIRST_INFECTION, REINFECTION, RECOVER_FIRST_INFECTION };
    case IR_STATE:
      return { BIRTHDEATH, FIRST_INFECTION, REINFECTION, RECOVER_REINFECTION };
    break;
    case R_STATE:
      return { BIRTHDEATH, FIRST_INFECTION, WANING };
    break;
    case P_STATE:
      return { BIRTHDEATH, FIRST_INFECTION, REINFECTION };
    break;
  }
}

std::vector<EventType> rate_changes(EventType ev) {
  switch(ev) {
    case FIRST_INFECTION:
      return { FIRST_INFECTION, REINFECTION, RECOVER_FIRST_INFECTION, BIRTHDEATH };
    case REINFECTION:
      return { FIRST_INFECTION, REINFECTION, RECOVER_REINFECTION };
    case RECOVER_FIRST_INFECTION:
      return { FIRST_INFECTION, REINFECTION, RECOVER_FIRST_INFECTION, WANING };
    case RECOVER_REINFECTION:
      return { FIRST_INFECTION, REINFECTION, RECOVER_REINFECTION, WANING };
    case WANING:
      return { REINFECTION, WANING };
  }
}

void update(std::vector<double> &states, Params p, EventType ev, std::vector<double> &rates, gsl_rng *rng){
   // which rates will need updates
  StateType from = ev != BIRTHDEATH ? 
    from(ev) :
    stateref[1+discrete_ran_weighted(rng, &states[1], NUM_OF_STATE_TYPES-1)];
  std::vector<EventType> updates = ev != BIRTHDEATH ?
    rate_changes(ev) :
    birthdeath_rate_changes(from);
      
  fromto(from(ev), to(ev), states);
  update_rates(states, p, rates, updates);
}