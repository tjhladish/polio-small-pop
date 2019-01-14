#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include <array>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <sstream>
#include <map>
#include <sys/time.h>

using namespace std;

#include "States.h"
#include "Params.h"
#include "Gillespie.h"
#include "Utils.h"
#include "PolioEvents.h"

enum ObsEventType {
  DETECTION, EXTINCTION, EVENT_MAX, TIME_MAX
};

const vector<const string> obsevent = {"D", "E", "V", "T"};

struct ObsEvent {
  unsigned int sim_id;
  double time;
  ObsEventType ot;
  ObsEvent(unsigned int si, double t, ObsEventType o) : sim_id(si), time(t), ot(o) {}
};

// TODO make these argument based
// string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
const string output_dir = "~/workspaces/polio-small-pop/output/";
//string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/Corrected SC Sims Results/";
const string ext = "_filter_0pcase_trunc.csv";

const string SEP = ","; // output separator--was ", "

const int SAMPLES    = 10; // Number of Simulations to run:
const int MAX_EVENTS = 1e8;
const int EPS_MAX    = 10;   // circulation interval considered for EPS calculation, in years

vector<int> allocate_data_structures(unsigned int len, int init=0) {
  return vector<int>(len, init);
}

int main(int argc, char* argv[]){
  string config = argc > 1 ? argv[1] : "main_gillespie.json";

  Params p = parseParams(config);
  //find expected compartment size
  const vector<double> compartments = equilibrium_fraction(&p);

  // simulation outcomes
  auto simevents = vector<ObsEvent>(); // used to create circulation intervals -- elements are actual time points
    
  //The Simulation
  for(int simid = 0; simid < SAMPLES; ++simid){
    //reset all parameters to original values after each run of the simulation
    double current_time = 0.0, last_time = 0.0;
    // kahan sum
    // double carry_time = 0.0, tmp_time = 0.0;

    auto rng = initialize_rng(simid);
    
    // initialize system
    vector<unsigned int> states;
    bool extinct;
    do { // require at least some infection to initialize sample
      states = multinomial_compartments(rng, compartments, p.Population[0]);
      extinct = (states[I1_STATE]+states[IR_STATE]) == 0;
    } while (extinct);
    
    auto rates = set_rates(states, &p);
    int nevents = 0;

    do {
      GEvent event = gillespie_ran_event(rng, rates);
      EventType ev = eventref[event.which];
      // this potentially has large + small machine precision error issues; kahan sum correction?
      // especial true larger time gets / smaller delta gets (i.e., larger populations)
      current_time = current_time + event.deltaT;
      
      // observation process
      if (ev == FIRST_INFECTION) { // && gillespie_ran_bool(rng, PIR*DET_RATE)
        // since these params subj of sensitivity study, but don't do anything to xmission
        // can analyze this in post
        simevents.push_back(ObsEvent(simid, current_time, DETECTION));
        last_time = current_time;
      }
      // states + rates
      update(states, &p, ev, rates, rng);
      
      extinct = (states[I1_STATE]+states[IR_STATE]) == 0;    
    } while (
      !extinct && // still possible to generate new cases
      nevents++ < MAX_EVENTS && // haven't exceeded event budget
      (current_time - last_time) < EPS_MAX // haven't gone too long without recordable event
    );
    
    // determine how simulation concluded
    ObsEventType last;
    if (extinct) {
      last = EXTINCTION;
    } else if (nevents >= MAX_EVENTS) {
      last = EVENT_MAX;
    } else if (current_time - last_time >= EPS_MAX) {
      last = TIME_MAX;
    } else {
      cerr << "SOMETHING AWRY " << simid << endl;
      exit(0); // some error condition?
    }
    simevents.push_back(ObsEvent(simid, current_time, last));
    
  }
  
  for (auto ev : simevents) {
    cout << ev.sim_id << ", " << ev.time << ", " << obsevent[ev.ot] << endl;
  }
  return 0;
}
