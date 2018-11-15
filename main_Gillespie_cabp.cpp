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

// string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
string output_dir = "~/workspaces/polio-small-pop/output/";
//string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/Corrected SC Sims Results/";
string ext = "_filter_0pcase_trunc.csv";

const string SEP = ","; // output separator--was ", "

#include "States.h"
#include "Params.h"
#include "Gillespie.h"
#include "PolioEvents.h"

enum OutputType {
  PCASE_TIME_OUT,
  PCASE_INCIDENCE_OUT,
  EXTINCTION_TIME_OUT,
  S_OUT,
  I1_OUT,
  R_OUT,
  P_OUT,
  IR_OUT,
  S_AT_PCASE_OUT,
  I1_AT_PCASE_OUT,
  R_AT_PCASE_OUT,
  P_AT_PCASE_OUT,
  IR_AT_PCASE_OUT,
  TIME_OUT,
  PCASE_TALLY_OUT,
  FIRST_INF_PER_YEAR_OUT,
  FIRST_INF_EVENT_TIMES_OUT,
  CIRCULATION_INTERVAL_OUT,
  NUM_OF_OUTPUT_TYPES
};

void output_results(vector<stringstream> &output_streams) {

    string base_filename = "stuff" + ext;
    vector<string> output_filenames(NUM_OF_OUTPUT_TYPES);

    output_filenames[PCASE_TIME_OUT           ] = output_dir + "time_of_pcases_N_"+ base_filename;
    output_filenames[PCASE_INCIDENCE_OUT      ] = output_dir + "num_p_cases_N_"+ base_filename;
    output_filenames[EXTINCTION_TIME_OUT      ] = output_dir + "TTE_N_"+ base_filename;
    output_filenames[S_OUT                    ] = output_dir + "S_"+ base_filename;
    output_filenames[I1_OUT                   ] = output_dir + "I1_"+ base_filename;
    output_filenames[R_OUT                    ] = output_dir + "R_"+ base_filename;
    output_filenames[P_OUT                    ] = output_dir + "P_"+ base_filename;
    output_filenames[IR_OUT                   ] = output_dir + "Ir_"+ base_filename;
    output_filenames[S_AT_PCASE_OUT           ] = output_dir + "S_at_pcase_"+ base_filename;
    output_filenames[I1_AT_PCASE_OUT          ] = output_dir + "I1_at_pcase_"+ base_filename;
    output_filenames[R_AT_PCASE_OUT           ] = output_dir + "R_at_pcase_"+ base_filename;
    output_filenames[P_AT_PCASE_OUT           ] = output_dir + "P_at_pcase_"+ base_filename;
    output_filenames[IR_AT_PCASE_OUT          ] = output_dir + "Ir_at_pcase_"+ base_filename;
    output_filenames[TIME_OUT                 ] = output_dir + "time_"+ base_filename;
    output_filenames[FIRST_INF_EVENT_TIMES_OUT] = output_dir + "first_inf_event_times_"+ base_filename;
    output_filenames[FIRST_INF_PER_YEAR_OUT   ] = output_dir + "first_inf_per_year_" + base_filename;
    output_filenames[CIRCULATION_INTERVAL_OUT ] = output_dir + "circulation_interval_"+base_filename;
    output_filenames[PCASE_TALLY_OUT          ] = output_dir + "pCases_per_year_" + base_filename;

//    for (int ot_idx = 0; ot_idx < NUM_OF_OUTPUT_TYPES; ++ot_idx) {
//        const OutputType ot = (OutputType) ot_idx;
        const OutputType ot = CIRCULATION_INTERVAL_OUT;
        ofstream ofs;
        ofs.open(output_filenames[ot]);
        ofs << output_streams[ot].rdbuf();
        ofs.close();
//    }
}

gsl_rng* initialize_rng(int seed = 0) {
  gsl_rng_env_setup(); // sets default rng from environment
  const gsl_rng_type* rng_type;
  rng_type = gsl_rng_default;
  gsl_rng* rng;
  rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(rng, 0);
  return rng;
}

map<OutputType, stringstream> init_output_streams(vector<OutputType> &outputs) {
  map<OutputType, stringstream> res();
  for(auto ot: outputs) res.put(ot, stringstream());
  return res;
}

enum ObsEventType {
  DETECTION, EXTINCTION, EVENT_MAX, TIME_MAX
}

struct ObsEvent {
  ObsEventType ot; double time;
}

int main(int argc, char* argv[]){
  string config = argc > 0 ? argv[0] : "main_gillespie.json";
  Params p = parseParams(config);

  map<OutputType, stringstream> output_streams = init_output_streams({ CIRCULATION_INTERVAL_OUT, PCASE_TIME_OUT });

  // TODO make these argument based
  const int SAMPLES    = 1000; // Number of Simulations to run:
  const int MAX_EVENTS = 1e8;
  const int EPS_RES    = 365;  // resolution of endemic potential statistic, in divisions per year
  const int EPS_MAX    = 20;   // circulation interval considered for EPS calculation, in years

  // simulation outcomes
  vector<double> case_detection_times; // used to create circulation intervals -- elements are actual time points
  // TODO output actual intervals - important once EPS_RES, EPS_MAX are program inputs?
  
  // record of observed intervals
  vector<int> extinction_intervals(EPS_MAX*EPS_RES, 0);
  vector<int> intercase_intervals(EPS_MAX*EPS_RES, 0);
  
  // record of intervals under case-at-t-0 assumption
  vector<int> initial_intercase_ints(EPS_MAX*EPS_RES, 0);
  vector<int> initial_extinction_ints(EPS_MAX*EPS_RES, 0);
  
  // censored intervals - i.e., unknown if extinction or case next
  vector<int> max_events_ints(EPS_MAX*EPS_RES, 0);
  vector<int> max_events_ints_no_icase(EPS_MAX*EPS_RES, 0);
  int lost_interest_no_icase = 0,
      lost_interest = 0;
    
  //find expected compartment size
  const vector<double> compartments = equilibrium_fraction(p);

  //The Simulation
  for(int i = 0; i < SAMPLES; ++i){
    //reset all parameters to original values after each run of the simulation
    double current_time = 0, last_time = 0;
    case_detection_times.clear();
    auto rng = initialize_rng(i);
    
    // initialize system
    std::vector<unsigned int> states;
    bool extinct;
    do { // require at least some infection to initialize sample
      states = multinomial_compartments(rng, compartments, p.Population);
      extinct = (states[I1_STATE]+states[IR_STATE]) == 0;
    } while (extinct);
    
    auto rates = set_rates(states, &p);
    int nevents = 0;

    do {
      GEvent event = gillespie_ran_event(rng, rates);
      EventType ev = eventref[event.which];
      // this potentially has large + small machine precision error issues; kahan sum correction?
      // especial true larger time gets / smaller delta gets (i.e., larger populations)
      current_time += event.deltaT;
      
      // observation process
      if (ev == FIRST_INFECTION && gillespie_ran_bool(rng, PIR*DET_RATE)) {
        case_detection_times.push_back(current_time);
        last_time = current_time;
      }
      update(&states, p, ev, &rates, rng);
      extinct = (states[I1_STATE]+states[IR_STATE]) == 0;    
    } while (
      !extinct && // still possible to generate new cases
      nevents++ < MAX_EVENTS && // haven't exceeded event budget
      (current_time - last_time) < EPS_MAX // haven't gone too long without recordable event
    );
    
    if (case_detection_times.size() > 0) { // saw at least one case
      // calculate intercase intervals
      // calculate initial interval
      
    } else { // saw no cases
      
    }
    
    // possible outcomes:
    // saw a case or not
    // if extinct:
    //  saw at least one case vs not
    // else (not extinct)
    //  ran out of events vs
    //  lost interest
    
    // intercase intervals == diffs of case detection times
    // extinction intervals == for sims ending in extinction, current_time - case_detection_times.back()
    // circulation intervals == sum of those two.
    // circulation can be done easily in post processing

      //run the simulation for 100 mil steps
      for(int j=0; j<MAX_EVENTS; ++j){

          //stopping condition
          if( (states[I1_STATE]+states[IR_STATE]) == 0.0 || time >= EPS_MAX){
              circInt.push_back(time); // not a good variable name--these aren't circulation intervals, they're case times
              for(unsigned int i = 0; i < circInt.size(); i++){
                  const double ci = i > 0 ? circInt[i] - circInt[i-1] : circInt[i];
                  const int eps_idx = ci < EPS_MAX ? (int) (ci*EPS_RES) : EPS_MAX*EPS_RES - 1;
                  eps_circ_ivls[eps_idx]++;
                  output_streams[CIRCULATION_INTERVAL_OUT] << circInt[i];
                  if(i < circInt.size() - 1){
                      eps_intercase_ivls[eps_idx]++;
                      output_streams[CIRCULATION_INTERVAL_OUT] << SEP;
                  }
              }
              output_streams[CIRCULATION_INTERVAL_OUT] << endl;
              
              for(unsigned int i = 0; i < pCaseDetection.size(); i++){
                  output_streams[PCASE_TIME_OUT] << pCaseDetection[i];
                  if(i < pCaseDetection.size() - 1){
                      output_streams[PCASE_TIME_OUT] << SEP;
                  }
              }
              output_streams[PCASE_TIME_OUT] << endl;
              break;
          }
      }

  }
  assert(eps_intercase_ivls.size() == eps_circ_ivls.size());
  for (unsigned int i = 0; i < eps_circ_ivls.size(); ++i) {
      cerr << eps_intercase_ivls[i] << SEP << eps_circ_ivls[i] << endl;
  }
  output_results(output_streams);
  return 0;
}
