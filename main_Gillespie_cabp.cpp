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

int main(){
    vector<stringstream> output_streams(NUM_OF_OUTPUT_TYPES);

    const int numSims=1000;                  // Number of Simulations to run:
    vector<double> TTE;                     // time to extinction vector -- use to get numerator for Eichner & Dietz stat
    vector<double> pCaseDetection;          // time of paralytic cases vector (also use for case-free periods)
    vector<double> totalParalyticCases;     // vector for num paralytic cases
    vector<double> pCasesPerYear;           // vector for paralytic cases per year
    vector<double> histogramCases(50,0);    // vector for counting number of cases per year
    vector<double> first_infections_per_year; //vector for calculating first infections per year
    vector<double> histogramFirstInf(1000,0); //vector for counting number of first infections per year
    vector<double> time_at_first_inf;       //time at first infection vector
    vector<double> circInt;                 //used to create circulation intervals -- elements are actual time points

    //int seed = 20;
    //mt19937 gen(seed);

    // random_device rd;                       // generates a random real number for the seed
    // mt19937 gen(rd());                      // random number generator

    auto rng = initialize_rng();

    const int MAX_EVENTS = 1e8;
    const int EPS_RES = 365; // resolution of endemic potential statistic, in divisions per year
    const int EPS_MAX = 20;  // circulation interval considered for EPS calculation, in years
    vector<int> eps_circ_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each
    vector<int> eps_intercase_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each

    //find expected compartment size
    const vector<double> compartments = equilibrium_fraction(p);

    //The Simulation
    for(int i=0; i<numSims; ++i){
        //reset all parameters to original values after each run of the simulation
        double time     = 0;
        double timeBet = 0;
        pCaseDetection.clear();
        pCasesPerYear.clear();
        first_infections_per_year.clear();
        time_at_first_inf.clear();
        circInt.clear();
        circInt.push_back(0);

        auto states = multinomial_compartments(rng, compartments, TOT);
        auto rates = set_rates(states, &p);

        //run the simulation for 100 mil steps
        for(int j=0; j<MAX_EVENTS; ++j){
            GEvent event = gillespie_ran_event(rng, rates);
            EventType ev = eventref[event.which];
            time += event.deltaT;
            //Pick the event that is to occur based on generated number and rate
            if (ev == FIRST_INFECTION && gillespie_ran_bool(rng, PIR*DET_RATE)) {
              if(timeBet > 0) pCaseDetection.push_back(time - timeBet);
              timeBet = time;
              circInt.push_back(time);
            }
            update(&states, p, ev, &rates, rng);
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
                        output_streams[CIRCULATION_INTERVAL_OUT]<< SEP;
                    }
                }
                output_streams[CIRCULATION_INTERVAL_OUT] << endl;
                
                for(unsigned int i = 0; i < pCaseDetection.size();i++){
                    output_streams[PCASE_TIME_OUT]<<pCaseDetection[i];
                    if(i < pCaseDetection.size() - 1){
                        output_streams[PCASE_TIME_OUT]<<SEP;
                    }
                }
                output_streams[PCASE_TIME_OUT]<<endl;
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
