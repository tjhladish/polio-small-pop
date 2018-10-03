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
#include <algorithm>


using namespace std;

string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/Corrected SC Sims Results/";
string ext = "_movement_test_const_pop1.csv";
const string SEP = ","; // output separator--was ", "

int villageCounter = 0;
uniform_real_distribution<> unifdis(0.0, 1.0);

enum StateType {S_STATE,
    I1_STATE,
    R_STATE,
    P_STATE,
    IR_STATE,
    NUM_OF_STATE_TYPES};

enum EventType {
    FIRST_INFECTION_EVENT,
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
    MOVEMENT,
    NUM_OF_EVENT_TYPES};

enum OutputType{
    S_OUT,
    I1_OUT,
    R_OUT,
    P_OUT,
    IR_OUT,
    NON_EXTINCT_VILLAGES,
    CIRCULATION_INTERVAL_OUT,
    NUM_OF_OUTPUT_TYPES
};


struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    vector<double> Tot;
};

struct Event{ //keeps track of which events are occurring to which village
    EventType event;
    int village;
};

const double KAPPA              = 0.4179; //waning depth parameter
const double RHO                = 0.2; //waning speed parameter

//other parameters
const int    numVillages        = 3;                   //total number of villages under consideration
//const vector<double> TOT        (numVillages,10000);  //total population size
const vector<double> TOT = {10000,5000,5000};
const double RECOVERY           = 13;    //recovery rate (individuals/year)
const double BETA               = 135;   //contact rate (individuals/year)
const double BIRTH              = 0.02; //birth rate (per year)
const double DEATH              = 0.02; //death rate (per year)
const double PIR                = 0.005;            //type 1 paralysis rate (naturally occurring cases)
const double DET_RATE           = 1.0;
const double MOVE_RATE          = 0.5; //needs to be updated at some point!!!!

vector<vector<double>> event_rates(numVillages, vector<double>(NUM_OF_EVENT_TYPES, 0.0));

vector<pair<double,double>> location;

int func_m(const gsl_vector * x, void * p, gsl_vector * f){
    Params * params = (Params *)p;
    const double recovery = (params->recovery);
    const double beta = (params->beta);
    const double birth = (params->birth);
    const double death = (params->death);
    const double kappa = (params->kappa);
    const double rho = (params->rho);
    const double Tot = (params->Tot[villageCounter]);
    const double S  = gsl_vector_get(x,S_STATE);
    const double I1 = gsl_vector_get(x,I1_STATE);
    const double R  = gsl_vector_get(x,R_STATE);
    const double P  = gsl_vector_get(x,P_STATE);
    const double Ir = gsl_vector_get(x,IR_STATE);
    
    gsl_vector_set (f, S_STATE,  birth*Tot - (beta*S*(I1+kappa*Ir))/Tot - death*S);
    gsl_vector_set (f, I1_STATE, (beta*S*(I1+kappa*Ir))/Tot-recovery*I1-death*I1);
    gsl_vector_set (f, R_STATE,  recovery*I1+(recovery/kappa)*Ir-rho*R-death*R);
    gsl_vector_set (f, P_STATE,  rho*R - (kappa*beta*P*(I1+kappa*Ir))/Tot - death*P);
    gsl_vector_set (f, IR_STATE, Tot - (S+I1+R+P+Ir));
    
    return GSL_SUCCESS;
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
    vector<double> compartments = {gsl_vector_get(workspace_F->x,S_STATE),
        gsl_vector_get(workspace_F->x,I1_STATE),
        gsl_vector_get(workspace_F->x,R_STATE),
        gsl_vector_get(workspace_F->x,P_STATE),
        gsl_vector_get(workspace_F->x,IR_STATE)};
    
    gsl_multiroot_fsolver_free(workspace_F);
    
    gsl_vector_free(x);
    
    return compartments;
}

vector<int> multinomial_Compartments(int num_Compartments,const vector<double> expectedComp,int i){
    
    const gsl_rng_type* T;
    gsl_rng* r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    //sets seed by time of day
    struct timeval tv;
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    //mySeed = 20;
    gsl_rng_set(r,mySeed);
    unsigned int num_Trials = TOT[i];
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

void initialize_rates(const double S, const double I1, const double R, const double P, const double Ir, const int i) {
    event_rates[i][FIRST_INFECTION_EVENT]                = BETA*S/TOT[i]*(I1+ KAPPA*Ir); //first infection event
    event_rates[i][REINFECTION_EVENT]                    = KAPPA*BETA*P/TOT[i]*(I1+ KAPPA*Ir); //reinfection event
    event_rates[i][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1; //first infected revovery event
    event_rates[i][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir; //reinfected recovery event
    event_rates[i][WANING_EVENT]                         = RHO*R; //waning event
    event_rates[i][BIRTH_EVENT]                          = (S+I1+R+P+Ir<TOT[i]) ? BIRTH*(S+I1+R+P+Ir) : 0; //birth event
    event_rates[i][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R; //natural death of fully immune individual
    event_rates[i][DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P; //natural death of partially susceptible individual
    event_rates[i][DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S; //nautral death of naive susceptible
    event_rates[i][DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir; //natural death of reinfected individual
    event_rates[i][DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1; //natural death of first infected
    event_rates[i][MOVEMENT]                             = MOVE_RATE*TOT[i]; //rate of movement from village i
}

vector<Event> sample_event(mt19937& gen, double& totalRate, const vector<double> S, const vector<double> I1, const vector<double> R, const vector<double> P, const vector<double> Ir,const vector<int> prevVillage, const double time,vector<double> &rateVec) {
    totalRate = 0.0;
    vector<Event> eventOccurrence;
    Event sampleEvent;
    for(int i = 0; i < numVillages; i++){
        if(find(prevVillage.begin(),prevVillage.end(),i)!=prevVillage.end() or time==0){//recalculate rates if something changed or at beginning of new sim
            rateVec[i] = 0;
            for(unsigned int k = 0; k < event_rates[i].size();k++){
                rateVec[i] += event_rates[i][k];
            }
        }
    }
    for(unsigned int i = 0; i < rateVec.size();i++){
        totalRate += rateVec[i];
    }
    double ran=totalRate*unifdis(gen);
    
    bool endLoop = false;
    
    for (int event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        for(int i = 0; i < numVillages; i++){
            if (choose_event(ran, event_rates[i][event])) {
                sampleEvent.event = (EventType) event;
                sampleEvent.village = i;
                eventOccurrence.push_back(sampleEvent);
                endLoop = true;
                break;
            }
        }
        if(endLoop){
            break;
        }
    }

    //++event_tally[event_type];
    return eventOccurrence;
}

void chooseCompartment(int moveInd, vector<double> &S, vector<double> &I1, vector<double> &R, vector<double> &P, vector<double> &Ir, int fromVil, int toVil){
    switch (moveInd) {
        case 0:
            S[fromVil]= --S[fromVil];
            S[toVil] = ++S[toVil];
            //change fromvil rates
            event_rates[fromVil][FIRST_INFECTION_EVENT]                = BETA*S[fromVil]/(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil])*(I1[fromVil]+ KAPPA*Ir[fromVil]);
            event_rates[fromVil][BIRTH_EVENT]                          = (S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]<TOT[fromVil]) ? BIRTH*(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]) : 0;
            event_rates[fromVil][DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S[fromVil];
            //change tovil rates
            event_rates[toVil][FIRST_INFECTION_EVENT]                = BETA*S[toVil]/(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil])*(I1[toVil]+ KAPPA*Ir[toVil]);
            event_rates[toVil][BIRTH_EVENT]                          = (S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]<TOT[toVil]) ? BIRTH*(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]) : 0;
            event_rates[toVil][DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S[toVil];
            break;
        case 1:
            I1[fromVil] = --I1[fromVil];
            I1[toVil] = ++I1[toVil];
            //change fromvil rates
            event_rates[fromVil][FIRST_INFECTION_EVENT]                = BETA*S[fromVil]/(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil])*(I1[fromVil]+ KAPPA*Ir[fromVil]);
            event_rates[fromVil][REINFECTION_EVENT]                    = KAPPA*BETA*P[fromVil]/(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil])*(I1[fromVil]+ KAPPA*Ir[fromVil]);
            event_rates[fromVil][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1[fromVil];
            event_rates[fromVil][BIRTH_EVENT]                          = (S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]<TOT[fromVil]) ? BIRTH*(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]) : 0;
            event_rates[fromVil][DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1[fromVil];
            //change tovil rates
            event_rates[toVil][FIRST_INFECTION_EVENT]                = BETA*S[toVil]/(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil])*(I1[toVil]+ KAPPA*Ir[toVil]);
            event_rates[toVil][REINFECTION_EVENT]                    = KAPPA*BETA*P[toVil]/(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil])*(I1[toVil]+ KAPPA*Ir[toVil]);
            event_rates[toVil][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1[toVil];
            event_rates[toVil][BIRTH_EVENT]                          = (S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]<TOT[toVil]) ? BIRTH*(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]) : 0;
            event_rates[toVil][DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1[toVil];
            break;
        case 2:
            R[fromVil] = --R[fromVil];
            R[toVil] = ++R[toVil];
            //change fromvil rates
            event_rates[fromVil][WANING_EVENT]                         = RHO*R[fromVil];
            event_rates[fromVil][BIRTH_EVENT]                          = (S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]<TOT[fromVil]) ? BIRTH*(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]) : 0;
            event_rates[fromVil][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R[fromVil];
            //change tovil rates
            event_rates[toVil][WANING_EVENT]                         = RHO*R[toVil];
            event_rates[toVil][BIRTH_EVENT]                          = (S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]<TOT[toVil]) ? BIRTH*(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]) : 0;
            event_rates[toVil][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R[toVil];
            break;
        case 3:
            P[fromVil] = --P[fromVil];
            P[toVil] = ++P[toVil];
            //change fromvil rates
            event_rates[fromVil][REINFECTION_EVENT]                    = KAPPA*BETA*P[fromVil]/(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil])*(I1[fromVil]+ KAPPA*Ir[fromVil]);
            event_rates[fromVil][BIRTH_EVENT]                          = (S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]<TOT[fromVil]) ? BIRTH*(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]) : 0;
            event_rates[fromVil][DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P[fromVil];
            //change tovil rates
            event_rates[toVil][REINFECTION_EVENT]                    = KAPPA*BETA*P[toVil]/(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil])*(I1[toVil]+ KAPPA*Ir[toVil]);
            event_rates[toVil][BIRTH_EVENT]                          = (S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]<TOT[toVil]) ? BIRTH*(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]) : 0;
            event_rates[toVil][DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P[toVil];
            break;
        case 4:
            Ir[fromVil] = --Ir[fromVil];
            Ir[toVil] = ++Ir[toVil];
            //change fromvil rates
            event_rates[fromVil][FIRST_INFECTION_EVENT]                = BETA*S[fromVil]/(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil])*(I1[fromVil]+ KAPPA*Ir[fromVil]);
            event_rates[fromVil][REINFECTION_EVENT]                    = KAPPA*BETA*P[fromVil]/(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil])*(I1[fromVil]+ KAPPA*Ir[fromVil]);
            event_rates[fromVil][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir[fromVil];
            event_rates[fromVil][BIRTH_EVENT]                          = (S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]<TOT[fromVil]) ? BIRTH*(S[fromVil]+I1[fromVil]+R[fromVil]+P[fromVil]+Ir[fromVil]) : 0;
            event_rates[fromVil][DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir[fromVil];
            //change tovil rates
            event_rates[toVil][FIRST_INFECTION_EVENT]                = BETA*S[toVil]/(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil])*(I1[toVil]+ KAPPA*Ir[toVil]);
            event_rates[toVil][REINFECTION_EVENT]                    = KAPPA*BETA*P[toVil]/(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil])*(I1[toVil]+ KAPPA*Ir[toVil]);
            event_rates[toVil][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir[toVil];
            event_rates[toVil][BIRTH_EVENT]                          = (S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]<TOT[toVil]) ? BIRTH*(S[toVil]+I1[toVil]+R[toVil]+P[toVil]+Ir[toVil]) : 0;
            event_rates[toVil][DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir[toVil];
            break;
        default:
            break;
    }
}

inline void process_first_infection_event(vector<double> &S, vector<double> &I1, const vector<double> R, const vector<double> P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    S[j]  = --S[j];
    I1[j] = ++I1[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]                  +R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S[j];
    event_rates[j][DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1[j];
}

inline void process_reinfection_event(const vector<double> S, const vector<double> I1, const vector<double> R, vector<double> &P, vector<double> &Ir, const int chosenVillage) {
    const int j = chosenVillage;
    P[j]  = --P[j];
    Ir[j] = ++Ir[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P[j];
    event_rates[j][DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir[j];
}

inline void process_recovery_from_first_infection_event(const vector<double> S, vector<double> &I1, vector<double> &R, const vector<double> P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    I1[j] = --I1[j];
    R[j]  = ++R[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1[j];
    event_rates[j][WANING_EVENT]                         = RHO*R[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R[j];
    event_rates[j][DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1[j];
}

inline void process_recovery_from_reinfection_event(const vector<double> S, const vector<double> I1, vector<double> &R, const vector<double> P, vector<double> &Ir,const int chosenVillage) {
    const int j = chosenVillage;
    Ir[j] = --Ir[j];
    R[j]  = ++R[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir[j];
    event_rates[j][WANING_EVENT]                         = RHO*R[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R[j];
    event_rates[j][DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir[j];
}

inline void process_waning_event(const vector<double> S, const vector<double> I1, vector<double> &R, vector<double>&P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    R[j] = --R[j];
    P[j] = ++P[j];
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][WANING_EVENT]                         = RHO*R[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R[j];
    event_rates[j][DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P[j];
}

inline void process_birth_event(vector<double> &S, const vector<double> I1, const vector<double> R, const vector<double> P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    S[j] = ++S[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S[j];
}

inline void process_death_from_recovered_event(const vector<double> S, const vector<double> I1, vector<double> &R, const vector<double> P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    R[j] = --R[j];
    event_rates[j][WANING_EVENT]                         = RHO*R[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_RECOVERED_EVENT]           = DEATH*R[j];
}

inline void process_death_from_partial_susceptible_event(const vector<double> S, const vector<double> I1, const vector<double> R, vector<double> &P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    P[j] = --P[j];
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT] = DEATH*P[j];
}

inline void process_death_from_fully_susceptible_event(vector<double> &S, const vector<double> I1, const vector<double> R, const vector<double> P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    S[j] = --S[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT]   = DEATH*S[j];
}

inline void process_death_from_reinfection_event(const vector<double> S, const vector<double> I1, const vector<double> R, const vector<double> P, vector<double> &Ir, const int chosenVillage) {
    const int j = chosenVillage;
    Ir[j] = --Ir[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_REINFECTION_EVENT]         = DEATH*Ir[j];
}

inline void process_death_from_first_infection_event(const vector<double> S, vector<double> &I1, const vector<double> R, const vector<double> P, const vector<double> Ir, const int chosenVillage) {
    const int j = chosenVillage;
    I1[j] = --I1[j];
    event_rates[j][FIRST_INFECTION_EVENT]                = BETA*S[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][REINFECTION_EVENT]                    = KAPPA*BETA*P[j]/TOT[j]*(I1[j]+ KAPPA*Ir[j]);
    event_rates[j][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1[j];
    event_rates[j][BIRTH_EVENT]                          = (S[j]+I1[j]+R[j]+P[j]+Ir[j]<TOT[j]) ? BIRTH*(S[j]+I1[j]+R[j]+P[j]+Ir[j]) : 0;
    event_rates[j][DEATH_FROM_FIRST_INFECTION_EVENT]     = DEATH*I1[j];
}

inline void process_movement_event(vector<double> &S, vector<double> &I1, vector<double> &R, vector<double> &P, vector<double> &Ir,
    const int chosenVillage, mt19937& gen, int toVil){
    const int fromVil = chosenVillage;
    double tempSum = 0;
    vector<double> tempVec = {S[fromVil],I1[fromVil],R[fromVil],P[fromVil],Ir[fromVil]};
    double tempTot = S[fromVil] + I1[fromVil] + R[fromVil] + P[fromVil] + Ir[fromVil];
    int moveInd = numeric_limits<int>::max();
    double rand = unifdis(gen);
    for(unsigned int k = 0; k < tempVec.size(); k++){
        tempSum+=tempVec[k];
        if(rand < (tempSum/tempTot)){
            moveInd = k;
            break;
        }
    }
    chooseCompartment(moveInd,S,I1,R,P,Ir,fromVil,toVil);
    double tempSum1 = 0;
    vector<double> tempVec1 = {S[toVil],I1[toVil],R[toVil],P[toVil],Ir[toVil]};
    double tempTot1 = S[toVil] + I1[toVil] + R[toVil] + P[toVil] + Ir[toVil];
    int moveInd1 = numeric_limits<int>::max();
    double rand1 = unifdis(gen);
    for(unsigned int k = 0; k < tempVec1.size(); k++){
        tempSum1+=tempVec1[k];
        if(rand1 < (tempSum1/tempTot1)){
            moveInd1 = k;
            break;
        }
    }
    const int newFrom = toVil;
    const int newTo = fromVil;
    chooseCompartment(moveInd1,S,I1,R,P,Ir,newFrom,newTo);
}

void output_results(vector<stringstream> &output_streams) {
    
    string numInEachVil = "";
    for(unsigned int i = 0; i < numVillages; i++){
        numInEachVil += to_string(int(TOT[i]));
    }
    
    string base_filename = numInEachVil+",beta_"+to_string(int(BETA))+",detect_rate_"+to_string(float(DET_RATE))+"rho_"+to_string(float(RHO))+ "numVillages_"+to_string(numVillages) + "migRate_"+to_string(float(MOVE_RATE)) + ext;
    //string base_filename = "test" + ext;
    vector<string> output_filenames(NUM_OF_OUTPUT_TYPES);
    
    output_filenames[NON_EXTINCT_VILLAGES     ] = output_dir + "non_extinct_villages_"+base_filename;
    output_filenames[CIRCULATION_INTERVAL_OUT ] = output_dir + "circulation_interval_"+base_filename;
    output_filenames[S_OUT                    ] = output_dir + "S_"+ base_filename;
    output_filenames[I1_OUT                   ] = output_dir + "I1_"+ base_filename;
    output_filenames[R_OUT                    ] = output_dir + "R_"+ base_filename;
    output_filenames[P_OUT                    ] = output_dir + "P_"+ base_filename;
    output_filenames[IR_OUT                   ] = output_dir + "Ir_"+ base_filename;

    
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
    
    //initialize size of vector of vector of compartments
    vector<vector<double>> compartments(numVillages);
    
    //initialize size of vector of vector of initial values
    vector<vector<int>> initialValues(numVillages);
    
    //initialize size of total rates for each village
    vector<double> rateVec(numVillages);
    
    //initialize size of vectors for individual compartments
    vector<double> S(numVillages);
    vector<double> I1(numVillages);
    vector<double> R(numVillages);
    vector<double> P(numVillages);
    vector<double> Ir(numVillages);
    
    
    vector<double> circInt; //used to create circulation intervals -- elements are actual time points
    
    vector<vector<double>> circulationInts(numVillages); //test vector for circulation intervals for each village
    
    vector<pair<double,int>> nonExtinctVillages; //(time,num villages)
    vector<bool> extinct(numVillages);
    int villageInCirc;
    
    uniform_int_distribution<> vilDis(0,(numVillages - 1));
    
    int seed = 20;
    mt19937 gen(seed);
    
    //random_device rd;                       // generates a random real number for the seed
    //mt19937 gen(rd());                      // random number generator
    
    const int numSims = 10000;

    //find expected compartment size for each village
    for(int i = 0; i < numVillages; i++){
        compartments[villageCounter] = initialize_compartments();
        villageCounter++;
    }

    //The Simulation
    for(int i = 0; i < numSims; i++){
        double time = 0;
        circInt.push_back(0);
        villageInCirc = numVillages;
        
        nonExtinctVillages.push_back(make_pair(0,villageInCirc));

        for(int i = 0; i < numVillages; i++){
            circulationInts[i].push_back(0);//first time for each vector for each village is 0
            extinct[i] = false;
        }
        
        vector<int> prevVillage{numeric_limits<int>::max(),numeric_limits<int>::max()};//initialize to max b/c village hasn't been chosen yet
        

        for(int i = 0; i < numVillages;i++){
            //set initial values for each village using multinomial dist
            initialValues[i] = multinomial_Compartments(compartments[i].size(),compartments[i],i);
            S[i]        = initialValues[i][S_STATE];   //naive susceptible (no previous contact w/virus, moves into I1)
            I1[i]       = initialValues[i][I1_STATE];  //first infected (only time paralytic case can occur, recovers into R)
            R[i]        = initialValues[i][R_STATE];   //recovered (fully immune, wanes into P)
            P[i]        = initialValues[i][P_STATE];   //partially susceptible (moves into Ir)
            Ir[i]       = initialValues[i][IR_STATE];  //reinfected (recovers into R)
            initialize_rates(S[i],I1[i],R[i],P[i],Ir[i],i);
        }
        for(int j = 0; j < 1e8; j++){
            double totalRate = 0;
            
            vector<Event> test = sample_event(gen, totalRate, S, I1, R, P, Ir, prevVillage, time, rateVec);
            assert(test.size()==1);
            EventType event_type = test[0].event;
            const int chosenVillage = test[0].village;
            prevVillage[0] = chosenVillage;//set equal for the next time an event is chosen
            //reset if movement isn't chosen
            if(event_type!=MOVEMENT){
                prevVillage[1] = numeric_limits<int>::max();
            }

            switch(event_type){
                case FIRST_INFECTION_EVENT: process_first_infection_event(S, I1, R, P, Ir, chosenVillage);
                {
                    double rr = unifdis(gen);
                    if(rr<(PIR*DET_RATE)){
                        circInt.push_back(time);
                        circulationInts[chosenVillage].push_back(time);
                        /*for(int i = 0; i < circulationInts[chosenVillage].size();i++){
                            cout<<"circulationInts["<<i<<"] "<<circulationInts[chosenVillage][i]<<"\n";
                        }
                        int nonChosenVillage = 0;
                        if(nonChosenVillage == chosenVillage){
                            nonChosenVillage = 1;
                        }
                        for(int j = 0; j < circulationInts[nonChosenVillage].size();j++){
                            cout<<"non chosen circulationInts["<<j<<"] "<<circulationInts[nonChosenVillage][j]<<"\n";
                        }*/
                    }
                }
                    break;
                case REINFECTION_EVENT:                     process_reinfection_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case RECOVERY_FROM_FIRST_INFECTION_EVENT:   process_recovery_from_first_infection_event(S, I1, R, P, Ir,chosenVillage);
                    break;
                case RECOVERY_FROM_REINFECTION_EVENT:       process_recovery_from_reinfection_event(S, I1, R, P, Ir,chosenVillage);
                    break;
                case WANING_EVENT:                          process_waning_event(S, I1, R, P, Ir, chosenVillage);                         break;
                case BIRTH_EVENT:                           process_birth_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case DEATH_FROM_RECOVERED_EVENT:            process_death_from_recovered_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT:  process_death_from_partial_susceptible_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT:    process_death_from_fully_susceptible_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case DEATH_FROM_REINFECTION_EVENT:          process_death_from_reinfection_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case DEATH_FROM_FIRST_INFECTION_EVENT:      process_death_from_first_infection_event(S, I1, R, P, Ir, chosenVillage);
                    break;
                case MOVEMENT:{
                    double tempSum = 0;
                    double rand1 = unifdis(gen);
                    double tempDenom = totalRate - rateVec[chosenVillage];
                    for(int k = 0; k < rateVec.size();k++){
                        if(k!=chosenVillage){
                            tempSum+=rateVec[k];
                            if(rand1 < (tempSum/(tempDenom))){
                                prevVillage[1] = k;
                                break;
                            }
                        }
                    }
                    process_movement_event(S,I1,R,P,Ir,chosenVillage,gen,prevVillage[1]);
                    break;
                }
                default:
                    cerr << "ERROR: Unsupported event type" << endl;
                    break;
            }
            for(int k = 0; k < numVillages; k++){
                assert(S[k] >= 0);
                assert(I1[k] >= 0);
                assert(R[k] >= 0);
                assert(P[k] >= 0);
                assert(Ir[k] >= 0);
            }
            //generate the time at which the event occurs
            exponential_distribution<>rng(totalRate);
            time+=rng(gen);
            for(int k = 0; k < numVillages; k++){
                if(I1[k]+Ir[k]==0 and extinct[k] == false){
                    villageInCirc--;
                    extinct[k] = true;
                    nonExtinctVillages.push_back(make_pair(time,villageInCirc));
                    break;
                }
            }
            bool zero_I1 = all_of(I1.begin(),I1.end(),[](int i){return i==0;});
            bool zero_Ir = all_of(Ir.begin(),Ir.end(),[](int i){return i==0;});
            

            //stopping condition
            if((zero_I1 and zero_Ir)){
                circInt.push_back(time);
                for(unsigned int i = 0; i < numVillages; i++){
                    circulationInts[i].push_back(time);
                }
                /*for(int i = 0; i < 2; i++){
                    for(int j = 0; j < circulationInts[i].size();j++){
                        cout<<"circulation ints ["<<i<<j<<"] "<<circulationInts[i][j]<<"\n";
                    }
                }*/
                for(unsigned int i = 0; i < nonExtinctVillages.size(); i++){
                    output_streams[NON_EXTINCT_VILLAGES] << nonExtinctVillages[i].first;
                    output_streams[NON_EXTINCT_VILLAGES] << SEP;
                    output_streams[NON_EXTINCT_VILLAGES] << nonExtinctVillages[i].second;
                    output_streams[NON_EXTINCT_VILLAGES] << endl;
                }
                for(unsigned int i = 0; i < numVillages;i++){
                    //cout<<"circulationInts["<<i<<"] size "<<circulationInts[i].size()<<"\n";
                    for(unsigned int j = 0; j < circulationInts[i].size();j++){
                        //cout<<"circulationInt["<<i<<"]["<<j<<"] "<<circulationInts[i][j]<<"\n";
                        output_streams[CIRCULATION_INTERVAL_OUT]<<circulationInts[i][j];
                        if(j < circulationInts[i].size() - 1){
                            output_streams[CIRCULATION_INTERVAL_OUT]<< SEP;
                        }
                    }
                    output_streams[CIRCULATION_INTERVAL_OUT] << endl;
                }
                for(unsigned int i = 0; i < numVillages; i++){
                    output_streams[S_OUT] << S[i] << SEP;
                    output_streams[I1_OUT] << I1[i] << SEP;
                    output_streams[R_OUT] << R[i] << SEP;
                    output_streams[P_OUT] << P[i] << SEP;
                    output_streams[IR_OUT] << Ir[i] << SEP;
                }
                output_streams[S_OUT] << endl;
                output_streams[I1_OUT] << endl;
                output_streams[R_OUT] << endl;
                output_streams[P_OUT] << endl;
                output_streams[IR_OUT] << endl;
                //output_streams[CIRCULATION_INTERVAL_OUT] << endl;
                /*for(unsigned int i = 0; i < circInt.size(); i++){
                    output_streams[CIRCULATION_INTERVAL_OUT]<<circInt[i];
                    if(i < circInt.size() - 1){
                        output_streams[CIRCULATION_INTERVAL_OUT]<< SEP;
                    }
                }*/
                //clear vector when done
                for(int i = 0; i < numVillages; i++){
                    for(unsigned int j = 0; j < circulationInts[i].size();j++){
                        circulationInts[i].clear();
                    }
                }
                nonExtinctVillages.clear();
                break;
            }
        }
    }
    output_results(output_streams);
    return 0;
}
