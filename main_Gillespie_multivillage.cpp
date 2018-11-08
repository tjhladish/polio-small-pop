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

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/Corrected SC Sims Results/";
string ext = "_test_EE1.csv";
const string SEP = ","; // output separator--was ", "

uniform_real_distribution<> unifdis(0.0, 1.0);

enum StateType {S_STATE,
                I1_STATE,
                R_STATE,
                P_STATE,
                IR_STATE,
                NUM_OF_STATE_TYPES};

enum EventType {FIRST_INFECTION_EVENT,
                REINFECTION_EVENT,
                RECOVERY_FROM_FIRST_INFECTION_EVENT,
                RECOVERY_FROM_REINFECTION_EVENT,
                WANING_EVENT,
                DEATH_EVENT,
                MOVE_EVENT,
                NUM_OF_EVENT_TYPES};

enum OutputType{//S_OUT,
                //I1_OUT,
                //R_OUT,
                //P_OUT,
                //IR_OUT,
                //NON_EXTINCT_VILLAGES,
                CIRCULATION_INTERVAL_OUT,
                EPIDEMIC_CURVE_OUT,
                NUM_OF_OUTPUT_TYPES};


struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    vector<double> Population;
};

struct VillageEvent{ //keeps track of which events are occurring to which village
    EventType event;
    int village;
};

//fast waning parameters:
//kappa = 0.4179
//rho = 0.2

//intermediate waning parameters:
//kappa = 0.6383
//rho = 0.04

//slow waning parameters:
//kappa = 0.8434
//rho = 0.02

const double KAPPA                 = 0.4179; //waning depth parameter
const double RHO                   = 0.2; //waning speed parameter

//other parameters
const vector<double> village_pop   = {10000};
const int numVillages              = village_pop.size(); //total number of villages under consideration
const int numDaysToRecover         = 28;
const double RECOVERY              = 365/numDaysToRecover;    //recovery rate (/year)
const double BETA                  = 135;   //contact rate (individuals/year)
const double lifespan              = 50;
const double BIRTH                 = 1/lifespan; //birth rate (per year)
const double DEATH                 = 1/lifespan; //death rate (per year)
const double PIR                   = 0.005;            //type 1 paralysis rate (naturally occurring cases)
const double DET_RATE              = 1.0;
const double expectedTimeUntilMove = 0; //years
const double MOVE_RATE             = expectedTimeUntilMove > 0 ? 1/expectedTimeUntilMove : 0;

vector<vector<double>> event_rates(numVillages, vector<double>(NUM_OF_EVENT_TYPES, 0.0));

// solve for endemic equilibrium of corresponding ODE system
int func_m(const gsl_vector * x, void * p, gsl_vector * f){
    Params * params = (Params *)p;
    const double recovery = (params->recovery);
    const double beta = (params->beta);
    const double birth = (params->birth);
    const double death = (params->death);
    const double kappa = (params->kappa);
    const double rho = (params->rho);
    const double Population = (params->Population[0]);
    const double S  = Population*gsl_vector_get(x,S_STATE);
    const double I1 = Population*gsl_vector_get(x,I1_STATE);
    const double R  = Population*gsl_vector_get(x,R_STATE);
    const double P  = Population*gsl_vector_get(x,P_STATE);
    const double Ir = Population*gsl_vector_get(x,IR_STATE);

    const double birthdt = birth*Population;
    const double foi = ((beta*(I1+kappa*Ir))/Population);
    const double infection1 = S*foi;
    const double infection2 = P*foi*kappa;
    const double recovery1 = recovery*I1;
    const double recovery2 = (recovery/kappa)*Ir;
    const double waning = rho*R;

    gsl_vector_set (f, S_STATE,  (birthdt - infection1 - death*S)/Population);
    gsl_vector_set (f, I1_STATE, (infection1 - recovery1 - death*I1)/Population);
    gsl_vector_set (f, R_STATE,  (recovery1 + recovery2 - waning - death*R)/Population);
    gsl_vector_set (f, P_STATE,  (waning - infection2 - death*P)/Population);
    gsl_vector_set (f, IR_STATE, (infection2 - recovery2 - death*Ir)/Population);

    return GSL_SUCCESS;
}

vector<double> initialize_compartment(int villageId) {
    //initial population from equilibrium values
    Params params ={};
    params.recovery = RECOVERY;
    params.beta = BETA;
    params.birth = BIRTH;
    params.death = DEATH;
    params.kappa = KAPPA;
    params.rho = RHO;
    params.Population = {village_pop[villageId]};

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
        gsl_vector_set(x,i,.1);
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

vector<int> multinomial_Compartments(int num_Compartments,const vector<double> expectedComp,int i,uint seed){
    const gsl_rng_type* T;
    gsl_rng* r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    //sets seed by time of day
    //struct timeval tv;
    //gettimeofday(&tv,0);
    unsigned long mySeed = seed;//tv.tv_sec + tv.tv_usec;
    //mySeed = 20;
    gsl_rng_set(r,mySeed);
    unsigned int num_Trials = village_pop[i];
    double *p = new double[num_Compartments];
    for(unsigned int i = 0; i < expectedComp.size(); i++){
        p[i] = num_Trials*expectedComp[i];
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

unsigned int rand_nonuniform_uint(const vector<double> weights, mt19937& gen) {
    const double totalWeight = accumulate(weights.begin(), weights.end(), 0.0);
    double ran = totalWeight*unifdis(gen);

    for (unsigned int idx = 0; idx < weights.size(); ++idx) {
        if (choose_event(ran, weights[idx])) {
            return idx;
            break;
        }
    }
    return weights.size(); // indicates failure to choose
}

void calculate_rates(const vector<double> &S, const vector<double> &I1, const vector<double> &R, const vector<double> &P, const vector<double> &Ir, const int i) {
    event_rates[i][FIRST_INFECTION_EVENT]                = BETA*S[i]*(I1[i]+ KAPPA*Ir[i])/village_pop[i]; //first infection event
    event_rates[i][REINFECTION_EVENT]                    = KAPPA*BETA*P[i]*(I1[i]+ KAPPA*Ir[i])/village_pop[i]; //reinfection event
    event_rates[i][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1[i]; //first infected revovery event
    event_rates[i][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*Ir[i]; //reinfected recovery event
    event_rates[i][WANING_EVENT]                         = RHO*R[i]; //waning event
    event_rates[i][DEATH_EVENT]                          = DEATH*village_pop[i]; //natural death
    event_rates[i][MOVE_EVENT]                           = MOVE_RATE*village_pop[i]; //rate of movement from village i
}

VillageEvent sample_event(mt19937& gen, double& totalRate, const vector<double> &S, const vector<double> &I1, const vector<double> &R, const vector<double> &P, const vector<double> &Ir, const double time) {
    totalRate = 0.0;
    VillageEvent ve;
    for(int i = 0; i < numVillages; i++){
        for(unsigned int k = 0; k < event_rates[i].size();k++){
            totalRate += event_rates[i][k];
        }
    }

    double ran = totalRate*unifdis(gen);

    for (int event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        for(int i = 0; i < numVillages; i++){
            if (choose_event(ran, event_rates[i][event])) {
                ve.event = (EventType) event;
                ve.village = i;
                return ve;
            }
        }
    }
    cerr << "Error: No event sampled" << endl;
    exit(-100);
}

void update_compartments(vector<double> &S, vector<double> &I1, vector<double> &R, vector<double> &P, vector<double> &Ir, const uint A, const uint B, const StateType state_from_A, const StateType state_from_B) {
    switch (state_from_A) {
        case S_STATE:  --S[A];  ++S[B];  break;
        case I1_STATE: --I1[A]; ++I1[B]; break;
        case R_STATE:  --R[A];  ++R[B];  break;
        case P_STATE:  --P[A];  ++P[B];  break;
        case IR_STATE: --Ir[A]; ++Ir[B]; break;
        default: break;
    }
    switch (state_from_B) {
        case S_STATE:  --S[B];  ++S[A];  break;
        case I1_STATE: --I1[B]; ++I1[A]; break;
        case R_STATE:  --R[B];  ++R[A];  break;
        case P_STATE:  --P[B];  ++P[A];  break;
        case IR_STATE: --Ir[B]; ++Ir[A]; break;
        default: break;
    }
    calculate_rates(S, I1, R, P, Ir, A);
    calculate_rates(S, I1, R, P, Ir, B);
}

inline void process_death_event(vector<double> &S, vector<double> &I1, vector<double> &R, vector<double> &P, vector<double> &Ir, const int chosenVillage, mt19937& gen) {
    const int j = chosenVillage;
    vector<double> rates(NUM_OF_STATE_TYPES, 0.0);
    rates[S_STATE]  = S[j];
    rates[I1_STATE] = I1[j];
    rates[R_STATE]  = R[j];
    rates[P_STATE]  = P[j];
    rates[IR_STATE] = Ir[j];

    StateType source_state = (StateType) rand_nonuniform_uint(rates, gen);

    if (source_state == S_STATE) return; // important to bail now, since nothing happens in this case

    ++S[j];
    switch(source_state) {
        case I1_STATE:
            --I1[j];
            break;
        case R_STATE:
            --R[j];
            break;
        case P_STATE:
            --P[j];
            break;
        case IR_STATE:
            --Ir[j];
            break;
        default:
            break;
    }
}

inline void process_movement_event(vector<double> &S, vector<double> &I1, vector<double> &R, vector<double> &P, vector<double> &Ir, const int A, const int B, mt19937& gen){
    vector<double> weights(NUM_OF_STATE_TYPES, 0.0);

    // Sample state of person to move from A to B
    weights[S_STATE]  = S[A];
    weights[I1_STATE] = I1[A];
    weights[R_STATE]  = R[A];
    weights[P_STATE]  = P[A];
    weights[IR_STATE] = Ir[A];

    StateType state_from_A = (StateType) rand_nonuniform_uint(weights, gen);

    // Sample state of person to move from B to A
    weights[S_STATE]  = S[B];
    weights[I1_STATE] = I1[B];
    weights[R_STATE]  = R[B];
    weights[P_STATE]  = P[B];
    weights[IR_STATE] = Ir[B];

    StateType state_from_B = (StateType) rand_nonuniform_uint(weights, gen);
    if (state_from_A == state_from_B) {
        return; // nothing actually happens
    } else {
        update_compartments(S, I1, R, P, Ir, A, B, state_from_A, state_from_B);
    }
}

void output_results(vector<stringstream> &output_streams) {
    string numInEachVil = "";
    for(int i = 0; i < numVillages; i++){
        numInEachVil += to_string(int(village_pop[i]));
    }

    string base_filename = numInEachVil+",beta_"+to_string(int(BETA))+",detect_rate_"+to_string(float(DET_RATE))+"rho_"+to_string(float(RHO))+ "numVillages_"+to_string(numVillages) + "migRate_"+to_string(float(MOVE_RATE)) + ext;
    vector<string> output_filenames(NUM_OF_OUTPUT_TYPES);

    output_filenames[CIRCULATION_INTERVAL_OUT ] = output_dir + "circulation_interval_"+base_filename;
    //output_filenames[EPIDEMIC_CURVE_OUT ] = output_dir + "epi_curve_"+base_filename;

    //for (int ot_idx = 0; ot_idx < NUM_OF_OUTPUT_TYPES; ++ot_idx) {
        //const OutputType ot = (OutputType) ot_idx;
        const OutputType ot = CIRCULATION_INTERVAL_OUT;
        ofstream ofs;
        ofs.open(output_filenames[ot]);
        ofs << output_streams[ot].rdbuf();
        ofs.close();
    //}
}

int main(){
    vector<stringstream> output_streams(NUM_OF_OUTPUT_TYPES);

    //initialize size of vector of vector of compartments
    vector<vector<double>> compartments(numVillages);

    //initialize size of vector of vector of initial values
    vector<vector<int>> initialValues(numVillages);

    //initialize size of total rates for each village
    vector<double> rateVec(numVillages);
    
    //count number of new infections (I1 infections) at time t to look at epidemic curve
    vector<pair<double,double>> epiCurve;

    //initialize size of vectors for individual compartments
    vector<double> S(numVillages, 0.0);
    vector<double> I1(numVillages, 0.0);
    vector<double> R(numVillages, 0.0);
    vector<double> P(numVillages, 0.0);
    vector<double> Ir(numVillages, 0.0);

    vector<double> circInt; //used to create circulation intervals -- elements are actual time points

    vector<vector<double>> circulationInts(numVillages); //test vector for circulation intervals for each village

    //vector<pair<double,int>> nonExtinctVillages; //(time,num villages)
    vector<bool> extinct(numVillages);
    int villageInCirc;

    //int seed = 20;
    //mt19937 gen(seed);

    random_device rd;                       // generates a random real number for the seed
    unsigned long int seed = rd();
    //uint seed = 2186064846;
    cerr << "seed: " << seed << endl;
    mt19937 gen(seed);                      // random number generator
    //mt19937 gen(rd());                      // random number generator

    const int numSims = 10000;

    //find expected compartment size for each village
    for(int i = 0; i < numVillages; i++){
        compartments[i] = initialize_compartment(i);
    }
    const int EPS_RES = 100; // resolution of endemic potential statistic, in divisions per year
    const int EPS_MAX = 50;  // circulation interval considered for EPS calculation, in years
    vector<int> eps_circ_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each
    vector<int> eps_intercase_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each

    //The Simulation
    for(int i = 0; i < numSims; i++){
        double time = 0;
        circInt.push_back(0);
        villageInCirc = numVillages;
        //epiCurve.clear();
        //nonExtinctVillages.push_back(make_pair(0,villageInCirc));

       for(int i = 0; i < numVillages; i++){
            circulationInts[i].push_back(0);//first time for each vector for each village is 0
            extinct[i] = false;
        }
         //cout<<"epi curve size after clear "<<epiCurve.size()<<"\n";
        //epiCurve.push_back(make_pair(0.0,1.0));//initialize at 1 new I1 infection
        for(int i = 0; i < numVillages;i++){
            //set initial values for each village using multinomial dist
            
            initialValues[i] = multinomial_Compartments(compartments[i].size(),compartments[i],i,gen());
            S[i]        = initialValues[i][S_STATE];   //naive susceptible (no previous contact w/virus, moves into I1)
            I1[i]       = initialValues[i][I1_STATE];  //first infected (only time paralytic case can occur, recovers into R)
            R[i]        = initialValues[i][R_STATE];   //recovered (fully immune, wanes into P)
            P[i]        = initialValues[i][P_STATE];   //partially susceptible (moves into Ir)
            Ir[i]       = initialValues[i][IR_STATE];  //reinfected (recovers into R)
            calculate_rates(S, I1, R, P, Ir, i);
        }

        for(int j = 0; j < 1e10; j++){
            double totalRate = 0;
            VillageEvent ve = sample_event(gen, totalRate, S, I1, R, P, Ir, time);
            EventType event_type = ve.event;
            const uint A = ve.village;

            switch(event_type){
                case FIRST_INFECTION_EVENT:                 --S[A]; ++I1[A];
                {
                    double rr = unifdis(gen);
                    if(rr<(PIR*DET_RATE)){
                        circInt.push_back(time);
                        circulationInts[A].push_back(time);
                    }
                }
                    break;
                case REINFECTION_EVENT:                     --P[A]; ++Ir[A];
                    break;
                case RECOVERY_FROM_FIRST_INFECTION_EVENT:   --I1[A]; ++R[A];
                    break;
                case RECOVERY_FROM_REINFECTION_EVENT:       --Ir[A]; ++R[A];
                    break;
                case WANING_EVENT:                          --R[A]; ++P[A];
                    break;
                case DEATH_EVENT:                           process_death_event(S, I1, R, P, Ir, A, gen);
                    break;
                case MOVE_EVENT:{
                    uint B = rand_nonuniform_uint(village_pop, gen);
                    process_movement_event(S,I1,R,P,Ir,A,B,gen);
                    break;
                }
                default:
                    cerr << "ERROR: Unsupported event type" << endl;
                    break;
            }
            calculate_rates(S,I1,R,P,Ir,A);

            for(int k = 0; k < numVillages; k++){
                assert(S[k] >= 0);
                assert(I1[k] >= 0);
                assert(R[k] >= 0);
                assert(P[k] >= 0);
                //if (Ir[k] < 0) cerr << S[k] << " " << I1[k] << " " << R[k] << " " << P[k] << " " << Ir[k] << " " << endl;
                assert(Ir[k] >= 0);
            }
            
            //generate the time at which the event occurs
            exponential_distribution<>rng(totalRate);
            time+=rng(gen);
            /*if(fmod(time,.01) < .001){
                epiCurve.push_back(make_pair(time, I1[A]));
            }*/
            /*for(int k = 0; k < numVillages; k++){
                if(I1[k]+Ir[k]==0 and extinct[k] == false){
                    villageInCirc--;
                    extinct[k] = true;
                    //nonExtinctVillages.push_back(make_pair(time,villageInCirc));
                    break;
                }
            }*/
            bool zero_I1 = all_of(I1.begin(),I1.end(),[](int i){return i==0;});
            bool zero_Ir = all_of(Ir.begin(),Ir.end(),[](int i){return i==0;});

            //stopping condition
            if((zero_I1 and zero_Ir)){
                //if(circulationInts[0].size() < 2){i--;}
                circInt.push_back(time);
                for(unsigned int i = 0; i < (unsigned) numVillages; i++){
                    circulationInts[i].push_back(time);
                }
                //cout<<"epiCurve size "<<epiCurve.size()<<"\n";
                /*for(unsigned int nI1Inf = 0; nI1Inf < (unsigned) epiCurve.size(); nI1Inf++){
                    cout<<"time "<<epiCurve[nI1Inf].first<<"\n";
                    cout<<"num inf "<<epiCurve[nI1Inf].second<<"\n";
                    output_streams[EPIDEMIC_CURVE_OUT]<<epiCurve[nI1Inf].first<<SEP<<epiCurve[nI1Inf].second<<endl;
                }*/
                
                for(unsigned int i = 0; i < (unsigned) numVillages;i++){
                    //cout<<"circulationInts["<<i<<"] size "<<circulationInts[i].size()<<"\n";
                    for(unsigned int j = 0; j < circulationInts[i].size();j++){
                        const double ci = j > 0 ? circulationInts[i][j] - circulationInts[i][j-1] : circulationInts[i][j];
                        const int eps_idx = ci < EPS_MAX ? (int) (ci*EPS_RES) : EPS_MAX*EPS_RES - 1;
                        eps_circ_ivls[eps_idx]++;
                        //cout<<"circulationInt["<<i<<"]["<<j<<"] "<<circulationInts[i][j]<<"\n";
                        output_streams[CIRCULATION_INTERVAL_OUT]<<circulationInts[i][j];
                        if(j < circulationInts[i].size() - 1){
                            eps_intercase_ivls[eps_idx]++;
                            output_streams[CIRCULATION_INTERVAL_OUT]<< SEP;
                        }
                    }
                    output_streams[CIRCULATION_INTERVAL_OUT] << endl;
                }
                //clear vector when done
                for(int i = 0; i < numVillages; i++){
                    circulationInts[i].clear();
                }
                break;
            }
        }
    }
    assert(eps_intercase_ivls.size() == eps_circ_ivls.size());
    for (unsigned int i = 0; i < eps_circ_ivls.size(); ++i) {
//        cerr << eps_intercase_ivls[i] << SEP << eps_circ_ivls[i] << endl;
    }
    output_results(output_streams);
    return 0;
}
