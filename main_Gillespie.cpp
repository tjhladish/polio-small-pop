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
#include<sys/time.h>

using namespace std;

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/Corrected SC Sims Results/";
string ext = "_ES_stat_multinomial.csv";

uniform_real_distribution<> unifdis(0.0, 1.0);

//fast waning parameters:
//kappa = 0.4179
//rho = 0.2

//intermediate waning parameters:
//kappa = 0.6383
//rho = 0.04

//slow waning parameters:
//kappa = 0.8434
//rho = 0.02

//waning parameters
const double KAPPA = 0.4179; //waning depth parameter
const double RHO = 0.2; //waning speed parameter

//other parameters
const double TOT        = 10000; //total population size
const double RECOVERY   = 13; //recovery rate (individuals/year)
const double BETA       = 135; //contact rate (individuals/year)
const double BIRTH      = 0.02; //birth rate (per year)
const double DEATH      = 0.02; //death rate (per year)
const double PIR        = 0.005; //type 1 paralysis rate (naturally occurring cases)
const double DET_RATE   = 1.0; //detection rate of paralytic case

enum EventType {FIRST_INFECTION_EVENT,
                REINFECTION_EVENT,
                RECOVERY_FROM_FIRST_INFECTION_EVENT,
                RECOVERY_FROM_REINFECTION_EVENT,
                WANING_EVENT,
                BIRTH_EVENT,
                DEATH_FROM_RECOVERED_EVENT,
                DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT,
                DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT,
                DEATH_FROM_REINFECTION_EVENT,
                DEATH_FROM_FIRST_INFECTION_EVENT,
                NUM_OF_EVENT_TYPES};

vector<double> event_rates(NUM_OF_EVENT_TYPES, 0.0);

enum OutputType {PCASE_INTERVAL_OUT,
                 PCASE_INCIDENCE_OUT,
                 EXTINCTION_TIME_OUT,
                 S_OUT,
                 I1_OUT,
                 R_OUT,
                 P_OUT,
                 IR_OUT,
                 TIME_OUT,
                 PCASE_TALLY_OUT,
                 FIRST_INFECTION_EVENT_TIME_OUT,
                 TIME_AT_FIRST_INF_OUT,
                 NUM_OF_OUTPUT_TYPES};

//const vector<double> kappas = {0.4179, 0.6383, 0.8434};//fast, intermed, slow
//const vector<double> rhos   = {0.2, 0.04, 0.02};//fast, intermed, slow

struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    double Tot;
};

//initializes population at equilibrium level of large population
int func_m(const gsl_vector * x, void * p, gsl_vector * f){
    Params * params = (Params *)p;
    const double recovery = (params->recovery);
    const double beta = (params->beta);
    const double birth = (params->birth);
    const double death = (params->death);
    const double kappa = (params->kappa);
    const double rho = (params->rho);
    const double Tot = (params->Tot);
    const double S = gsl_vector_get(x,0);
    const double I1 = gsl_vector_get(x,1);
    const double R = gsl_vector_get(x,2);
    const double P = gsl_vector_get(x,3);
    const double Ir = gsl_vector_get(x,4);


    gsl_vector_set (f,0,birth*Tot - (beta*S*(I1+kappa*Ir))/Tot - death*S);
    gsl_vector_set (f,1,(beta*S*(I1+kappa*Ir))/Tot-recovery*I1-death*I1);
    gsl_vector_set (f,2,recovery*I1+(recovery/kappa)*Ir-rho*R-death*R);
    gsl_vector_set (f,3,rho*R - (kappa*beta*P*(I1+kappa*Ir))/Tot - death*P);
    gsl_vector_set (f,4,Tot - (S+I1+R+P+Ir));

    return GSL_SUCCESS;
}

vector<int> multinomial_Compartments(int num_Compartments,const vector<double> expectedComp){
    
    const gsl_rng_type* T;
    gsl_rng* r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    //sets seed by time of day
    struct timeval tv;
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    gsl_rng_set(r,mySeed);
    unsigned int num_Trials = TOT;
    double *p = new double[num_Compartments];
    for(unsigned int i = 0; i < expectedComp.size(); i++){
        p[i] = expectedComp[i];
    }
    vector<int> initialCompartments(num_Compartments);
    //generates weights for compartments using equilibrium value from large population

    unsigned int *n = new unsigned int[num_Compartments];
    gsl_ran_multinomial(r, num_Compartments,num_Trials,p,n);
    for(int i = 0; i < num_Compartments; i++){
        initialCompartments[i] = n[i];
    }
    delete[] n;
    delete[] p;
    gsl_rng_free(r);
    return initialCompartments;
}

bool choose_event(double &ran, const double p) {
    if (ran < p) {
        return true;
    } else {
        ran -= p;
        return false;
    }
}

void initialize_rates (const double S, const double I1, const double R, const double P, const double Ir) {
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir); //first infection event
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir); //reinfection event
    event_rates[RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1; //first infected revovery event
    event_rates[RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir; //reinfected recovery event
    event_rates[WANING_EVENT]                         = RHO*R; //waning event
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0; //birth event
    event_rates[DEATH_FROM_RECOVERED_EVENT]           = DEATH*R; //natural death of fully immune individual
    event_rates[DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P; //natural death of partially susceptible individual
    event_rates[DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S; //nautral death of naive susceptible
    event_rates[DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir; //natural death of reinfected individual
    event_rates[DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1; //natural death of first infected
}


inline void process_first_infection_event                (double &S, double &I1, const double R, const double P, const double Ir) {
    --S; ++I1;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S;
    event_rates[DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1;
}

inline void process_reinfection_event                    (const double S, const double I1, const double R, double &P, double &Ir) {
    --P; ++Ir;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P;
    event_rates[DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir;
}

inline void process_recovery_from_first_infection_event  (const double S, double &I1, double &R, const double P, const double Ir) {
    --I1; ++R;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1;
    event_rates[WANING_EVENT]                         = RHO*R;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_RECOVERED_EVENT]           = DEATH*R;
    event_rates[DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1;
}

inline void process_recovery_from_reinfection_event      (const double S, const double I1, double &R, const double P, double &Ir) {
    --Ir; ++R;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir;
    event_rates[WANING_EVENT]                         = RHO*R;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_RECOVERED_EVENT]           = DEATH*R;
    event_rates[DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir;
}

inline void process_waning_event                         (const double S, const double I1, double &R, double &P, const double Ir) {
    --R;  ++P;
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[WANING_EVENT]                         = RHO*R;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_RECOVERED_EVENT]           = DEATH*R;
    event_rates[DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P;
}

inline void process_birth_event                          (double &S, const double I1, const double R, const double P, const double Ir) {
    ++S;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S;
}

inline void process_death_from_recovered_event           (const double S, const double I1, double &R, const double P, const double Ir) {
    --R;
    event_rates[WANING_EVENT]                         = RHO*R;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_RECOVERED_EVENT]           = DEATH*R;
}

inline void process_death_from_partial_susceptible_event (const double S, const double I1, const double R, double &P, const double Ir) {
    --P;
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P;
}

inline void process_death_from_fully_susceptible_event   (double &S, const double I1, const double R, const double P, const double Ir) {
    --S;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S;
}

inline void process_death_from_reinfection_event         (const double S, const double I1, const double R, const double P, double &Ir) {
    --Ir;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir;
}

inline void process_death_from_first_infection_event     (const double S, double &I1, const double R, const double P, const double Ir) {
    --I1;
    event_rates[FIRST_INFECTION_EVENT]                = BETA*S/TOT*(I1+ KAPPA*Ir);
    event_rates[REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT*(I1+ KAPPA*Ir);
    event_rates[RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1;
    event_rates[BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT) ? BIRTH*(S+I1+R+P+Ir) : 0;
    event_rates[DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1;
}


EventType sample_event(mt19937& gen, double& totalRate, const double S, const double I1, const double R, const double P, const double Ir) {
    totalRate = 0.0;
    for (auto rate: event_rates) totalRate += rate;

    //generate unifrn
    double ran=totalRate*unifdis(gen);

    for (int event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        if (choose_event(ran, event_rates[event])) { return (EventType) event; }
    }
    //++event_tally[event_type];*/
    return NUM_OF_EVENT_TYPES;
}

vector<double> initialize_compartments() {
    //initial population from equilibrium values
    Params params = {RECOVERY, BETA, BIRTH, DEATH, KAPPA, RHO, TOT};

    int i, times, status;
    gsl_multiroot_function F;
    gsl_multiroot_fsolver *workspace_F;
    gsl_vector *x;
    int num_dimensions = 5;

    x = gsl_vector_alloc(num_dimensions);

    workspace_F = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,num_dimensions);
    printf("F solver: %s\n", gsl_multiroot_fsolver_name(workspace_F));
    F.f=&func_m;
    F.n=num_dimensions;
    F.params = &params;
    int MAXTIMES = 100;
    /* set initial value */
    for(i = 0; i < num_dimensions; i++){
        gsl_vector_set(x,i,100);
    }

    /* set solver */
    gsl_multiroot_fsolver_set(workspace_F,&F, x);

    /* main loop */
    for(times = 0; times < MAXTIMES; times++) {
        status = gsl_multiroot_fsolver_iterate(workspace_F);

        /*fprintf(stderr, "%d times: ", times);
        for(i = 0; i < num_dimensions; i++) {
            fprintf(stderr, "%10.3e ", gsl_vector_get(workspace_F->x, i));
        }
        fprintf(stderr, "\n");
        */
        if((status == GSL_EBADFUNC) || (status == GSL_ENOPROG))
        {
            fprintf(stderr, "Status: %s\n", gsl_strerror(status));
            break;
        }
    }

    assert(status == GSL_ENOPROG);
    /* print answer */
    for(i = 0; i < num_dimensions; i++) {
        fprintf(stderr, "%3d %25.17e\n", i, gsl_vector_get(workspace_F->x, i));
    }

    assert(num_dimensions==5);
    vector<double> compartments =   {gsl_vector_get(workspace_F->x,0),
                                     gsl_vector_get(workspace_F->x,1),
                                     gsl_vector_get(workspace_F->x,2),
                                     gsl_vector_get(workspace_F->x,3),
                                     gsl_vector_get(workspace_F->x,4)};
    
    /*map<string, double> compartments = {{"S",  gsl_vector_get(workspace_F->x, 0)},
                                        {"I1", gsl_vector_get(workspace_F->x, 1)},
                                        {"R",  gsl_vector_get(workspace_F->x, 2)},
                                        {"P",  gsl_vector_get(workspace_F->x, 3)},
                                        {"Ir", gsl_vector_get(workspace_F->x, 4)}};*/

    gsl_multiroot_fsolver_free(workspace_F);

    /* free x */
    gsl_vector_free(x);

    return compartments;
}

void output_results(vector<stringstream> &output_streams) {

    string base_filename = to_string(TOT)+",beta_"+to_string(BETA)+",detect_rate_"+to_string(DET_RATE)+"rho_"+to_string(RHO)+ ext;
    map<OutputType,string> output_filenames = { {PCASE_INTERVAL_OUT,  output_dir + "time_between_pcases_N_"+ base_filename},
                                                {PCASE_INCIDENCE_OUT, output_dir + "num_p_cases_N_"+ base_filename        },
                                                {EXTINCTION_TIME_OUT, output_dir + "TTE_N_"+ base_filename                },
                                                {S_OUT,               output_dir + "S_"+ base_filename                    },
                                                {I1_OUT,              output_dir + "I1_"+ base_filename                   },
                                                {R_OUT,               output_dir + "R_"+ base_filename                    },
                                                {P_OUT,               output_dir + "P_"+ base_filename                    },
                                                {IR_OUT,              output_dir + "Ir_"+ base_filename                   },
                                                {TIME_OUT,            output_dir + "time_"+ base_filename                 },
                                                {FIRST_INFECTION_EVENT_TIME_OUT, output_dir + "first_inf_event_time_ "+ base_filename
                                                    },
                                                {TIME_AT_FIRST_INF_OUT, output_dir + "time_at_first_inf_" + base_filename},
                                                {PCASE_TALLY_OUT,     output_dir + "pCases_per_year_" + base_filename     }};


    for (int ot_idx = 0; ot_idx < NUM_OF_OUTPUT_TYPES; ++ot_idx) {
        const OutputType ot = (OutputType) ot_idx;
        ofstream ofs;
        ofs.open(output_filenames[ot]);
        ofs << output_streams[ot].rdbuf();
        ofs.close();
    }
}

int main(){
    vector<stringstream> output_streams(NUM_OF_OUTPUT_TYPES);

    const int numSims=10000;                  // Number of Simulations to run:
    vector<double> TTE;                     // time to extinction vector -- use to get numerator for Eichner & Dietz stat
    vector<double> pCaseDetection;          // time between paralytic cases vector (also use for case-free periods)
    vector<double> totalParalyticCases;     // vector for num paralytic cases
    vector<double> pCasesPerYear;           // vector for paralytic cases per year
    vector<double> histogramCases(50,0);    // vector for counting number of cases per year
    vector<double> first_infections_per_year; //vector for calculating first infections per year
    vector<double> histogramFirstInf(1000,0); //vector for counting number of first infections per year
    vector<double> time_at_first_inf;       //time at first infection vector

    //int seed = 20;
    //mt19937 gen(seed);

    random_device rd;                       // generates a random real number for the seed
    mt19937 gen(rd());                      // random number generator
    
    //find expected compartment size
    const vector<double> compartments = initialize_compartments();

    /*const int S_initial  = compartments.at("S"); //naive susceptible (no previous contact w/virus, moves into I1)
    const int I1_initial = compartments.at("I1"); //first infected (only time paralytic case can occur, recovers into R)
    const int R_initial  = compartments.at("R"); //recovered (fully immune, wanes into P)
    const int P_initial  = compartments.at("P"); //partially susceptible (moves into Ir)
    const int Ir_initial = compartments.at("Ir"); //reinfected (recovers into R)*/

    //The Simulation
    for(int i=0;i<numSims;++i){
        //reset all parameters to original values after each run of the simulation
        vector<int> initialValues = multinomial_Compartments(compartments.size(),compartments);
        double S        = initialValues[0];
        double I1       = initialValues[1];
        double R        = initialValues[2];
        double P        = initialValues[3];
        double Ir       = initialValues[4];
        double tsc      = 0; //used to calculate time between detected paralytic cases
        double time     = 0;
        int countPIR    = 0;
        pCaseDetection.clear();
        pCasesPerYear.clear();
        first_infections_per_year.clear();
        time_at_first_inf.clear();
        initialize_rates(S, I1, R, P, Ir);

        //run the simulation for 1 mill steps
        for(int j=0;j<1000001;++j){
            double totalRate = 0;
            EventType event_type = sample_event(gen, totalRate, S, I1, R, P, Ir);
            //Pick the event that is to occur based on generated number and rate
            switch (event_type) {
                case FIRST_INFECTION_EVENT: process_first_infection_event(S, I1, R, P, Ir);
                    //generate a unif random real num to determine if a paralytic case is detected
                    {
                        const int year = (int) time;
                        if((unsigned) year >= first_infections_per_year.size()){
                            first_infections_per_year.resize(year + 1,0);
                        }
                        first_infections_per_year[year]++;
                        time_at_first_inf.push_back(time);
                        
                        double rr = unifdis(gen);
                        if(rr<(PIR*DET_RATE)){
                            countPIR++;
                            //if(countPIR > 1){ //comment out to calculate Eichner & Dietz statistic
                                pCaseDetection.push_back(time-tsc);
                                //tsc=time;
                            //}
                            if((unsigned) year >= pCasesPerYear.size()){
                                pCasesPerYear.resize(year + 1, 0);
                            }
                            pCasesPerYear[year]++;
                            tsc=time;
                        }
                    }
                    break;
                case REINFECTION_EVENT:                     process_reinfection_event(S, I1, R, P, Ir);
                    break;
                case RECOVERY_FROM_FIRST_INFECTION_EVENT:   process_recovery_from_first_infection_event(S, I1, R, P, Ir);  break;
                case RECOVERY_FROM_REINFECTION_EVENT:       process_recovery_from_reinfection_event(S, I1, R, P, Ir);      break;
                case WANING_EVENT:                          process_waning_event(S, I1, R, P, Ir);                         break;
                case BIRTH_EVENT:                           process_birth_event(S, I1, R, P, Ir);                          break;
                case DEATH_FROM_RECOVERED_EVENT:            process_death_from_recovered_event(S, I1, R, P, Ir);           break;
                case DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT:  process_death_from_partial_susceptible_event(S, I1, R, P, Ir); break;
                case DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT:    process_death_from_fully_susceptible_event(S, I1, R, P, Ir);   break;
                case DEATH_FROM_REINFECTION_EVENT:          process_death_from_reinfection_event(S, I1, R, P, Ir);         break;
                case DEATH_FROM_FIRST_INFECTION_EVENT:      process_death_from_first_infection_event(S, I1, R, P, Ir);     break;
                default:
                    cerr << "ERROR: Unsupported event type" << endl;
                    break;
            }

            //generate the time at which the event occurs
            exponential_distribution<>rng(totalRate);
            time+=rng(gen);

            output_streams[S_OUT] << S      << ", ";
            output_streams[I1_OUT] << I1     << ", ";
            output_streams[R_OUT] << R      << ", ";
            output_streams[P_OUT] << P      << ", ";
            output_streams[IR_OUT] << Ir     << ", ";
            output_streams[TIME_OUT] << time    << ", ";

            //stopping condition
            if((I1+Ir==0) or time>10){
                totalParalyticCases.push_back(countPIR);
                TTE.push_back(time);
                pCaseDetection.push_back(time-tsc);//use for Eichner & Dietz statistic
                for(unsigned int i = 0; i < time_at_first_inf.size(); i++){
                    output_streams[TIME_AT_FIRST_INF_OUT]<<time_at_first_inf[i];
                    if(i < time_at_first_inf.size() - 1){
                        output_streams[TIME_AT_FIRST_INF_OUT]<< ", ";
                    }
                }
                output_streams[TIME_AT_FIRST_INF_OUT] <<"\n";
                const double fractional_year = time - (int) time;
                if (fractional_year > 0) {
                    pCasesPerYear.resize((int) time + 1, 0);
                    first_infections_per_year.resize((int) time + 1,0);
                }
                for (auto count: pCasesPerYear) histogramCases[count]++;
                for(unsigned int i = 0; i < first_infections_per_year.size();i++){
                    output_streams[FIRST_INFECTION_EVENT_TIME_OUT]<< first_infections_per_year[i];
                    if(i < first_infections_per_year.size()-1){
                        output_streams[FIRST_INFECTION_EVENT_TIME_OUT]<<",";
                    }
                }
                output_streams[FIRST_INFECTION_EVENT_TIME_OUT]<<"\n";
                /*for (auto count_i: first_infections_per_year){
                    assert(count_i < histogramFirstInf.size());
                    histogramFirstInf[count_i]++;
                }*/
                if (time != (int) time) {
                    const int last_years_count = pCasesPerYear.back();
                    
                    if (last_years_count != 0) {
                        histogramCases[last_years_count] -= 1.0 - fractional_year;
                    } else {
                        histogramCases[0] += fractional_year;
                    }
                    
                    /*const int last_years_count_inf = first_infections_per_year.back();
                    if(last_years_count_inf!=0){
                        histogramFirstInf[last_years_count_inf] -= 1.0 - fractional_year;
                    }
                    else{
                        histogramFirstInf[0] += fractional_year;
                    }*/
                }
                /*if(pCasesPerYear.size() == 1){
                    pCasesPerYear[0] = time;
                }*/
                //cout<<"pcases size "<<pCasesPerYear.size()<<endl;

                /*cout << "end time: " << time << "\npcases time series: ";
                for(unsigned int i = 0; i < pCasesPerYear.size() - 1; i++){
                    cout<<pCasesPerYear[i]<<",";
                }
                if (pCasesPerYear.size()>0) cout<<pCasesPerYear.back() << endl;
                cout << "pcase tally: ";
                */
                if (i == numSims-1) {
                    for(double ptally: histogramCases){ output_streams[PCASE_TALLY_OUT] << ptally << ",";} output_streams[PCASE_TALLY_OUT] << endl;
                    //for(double itally: histogramFirstInf){ output_streams[FIRST_INFECTION_EVENT_TIME_OUT] << itally << ",";}
                     //   output_streams[FIRST_INFECTION_EVENT_TIME_OUT] << endl;
                }
                //for(double pcase: pCaseDetection){output_streams[PCASE_INTERVAL_OUT]<<pcase<<",";}
                //output_streams[PCASE_INTERVAL_OUT]<<endl;
                for(unsigned int i = 0; i < pCaseDetection.size(); i++){
                    output_streams[PCASE_INTERVAL_OUT] << pCaseDetection[i];
                    if(i < pCaseDetection.size()-1){
                        output_streams[PCASE_INTERVAL_OUT]<<",";
                    }
                }
                output_streams[PCASE_INTERVAL_OUT] << "\n";
                output_streams[S_OUT] << S <<"\n";
                output_streams[I1_OUT] << I1 << " \n ";
                output_streams[R_OUT] << R << " \n ";
                output_streams[P_OUT] << P << " \n ";
                output_streams[IR_OUT] << Ir << " \n ";
                output_streams[TIME_OUT] << time << " \n";
                break;
            }
        }

    }
    for(unsigned int i = 0; i < totalParalyticCases.size(); i++){
        output_streams[PCASE_INCIDENCE_OUT] <<totalParalyticCases[i]<<"\n";
    }
    for (unsigned int i = 0; i < TTE.size(); i++) {
        output_streams[EXTINCTION_TIME_OUT] <<TTE[i]<<"\n";
    }
    output_results(output_streams);
    return 0;
}
