#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include <array>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
//#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <sstream>
#include <map>

using namespace std;

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/Corrected SC Sims Results/";
string ext = "_corrected_dist_1.csv";

uniform_real_distribution<> unifdis(0.0, 1.0);

//waning parameters
const double KAPPA = 0.4179; //kappas[atoi(argv[1])];
const double RHO = 0.2; //rhos[atoi(argv[1])];

//other parameters
const double TOT = 10000;
const double RECOVERY = 13;//gamma
const double BETA = 135;
const double BIRTH = 0.02;
const double DEATH = 0.02;
const double PIR =0.005; //type 1 paralysis rate
const double DET_RATE = 1.0;

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
                 NUM_OF_OUTPUT_TYPES };

const vector<double> kappas = {0.4179, 0.6383, 0.8434};//fast, intermed, slow
const vector<double> rhos   = {0.2, 0.04, 0.02};//fast, intermed, slow

struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    double Tot;
};

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
    //gsl_vector_set (f,4,kappa*beta*P*(I1+kappa*Ir)/Tot-(recovery/kappa)*Ir - death*Ir);
    
    
    return GSL_SUCCESS;
    
}

EventType sample_event(mt19937& gen, double& totalRate, const double S1, const double I11, const double R1, const double P1, const double Ir1) {
    //generate unifrn
    double ran=unifdis(gen);

    //update the transition rates
    double birthRate;
    if((S1+I11+R1+P1+Ir1<TOT)){
        birthRate = BIRTH*(S1+I11+R1+P1+Ir1);
    }
    else{
        birthRate = 0;
    }
    const double infect1 = BETA*S1/TOT*(I11+ KAPPA*Ir1);
    const double recover1 = RECOVERY*I11;
    const double wane = RHO*R1;
    const double infectr = KAPPA*BETA*P1/TOT*(I11+ KAPPA*Ir1);
    const double recover2 = (RECOVERY/KAPPA)*Ir1;
    const double deathS = DEATH*S1;
    const double deathI1 = DEATH*I11;
    const double deathR = DEATH*R1;
    const double deathP = DEATH*P1;
    const double deathIr = DEATH*Ir1;

    totalRate = birthRate+infect1+recover1+wane+infectr+recover2+deathS+deathI1+deathR+deathP+deathIr;

    EventType event_type;
    if(ran< infect1/totalRate){
        event_type = FIRST_INFECTION_EVENT;
    } else if(ran<((infectr+infect1)/totalRate)){
        event_type = REINFECTION_EVENT;
    } else if(ran<((recover1+infectr+infect1)/totalRate)){
        event_type = RECOVERY_FROM_FIRST_INFECTION_EVENT;
    } else if(ran<((recover2+recover1+infectr+infect1)/totalRate)){
        event_type = RECOVERY_FROM_REINFECTION_EVENT;
    } else if(ran<((wane+recover2+recover1+infectr+infect1)/totalRate)){
        event_type = WANING_EVENT;
    } else if(ran<((birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
        event_type = BIRTH_EVENT;
    } else if(ran<((deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
        event_type = DEATH_FROM_RECOVERED_EVENT;
    } else if(ran<((deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
        event_type = DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT;
    } else if(ran<((deathS+deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
        event_type = DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT;
    } else if(ran<((deathIr+deathS+deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
        event_type = DEATH_FROM_REINFECTION_EVENT;
    } else{
        event_type = DEATH_FROM_FIRST_INFECTION_EVENT;
    }
    return event_type;
}

map<string,double> initialize_compartments() {
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
    gsl_multiroot_fsolver_set(workspace_F, &F, x);
    
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
    map<string, double> compartments = {{"S",  gsl_vector_get(workspace_F->x, 0)},
                                        {"I1", gsl_vector_get(workspace_F->x, 1)},
                                        {"R",  gsl_vector_get(workspace_F->x, 2)},
                                        {"P",  gsl_vector_get(workspace_F->x, 3)},
                                        {"Ir", gsl_vector_get(workspace_F->x, 4)}};
   
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

    const int numSims=100;                  // Number of Simulations to run:
    vector<double> TTE;                     // time to extinction vector
    vector<double> pCaseDetection;          // time between paralytic cases vector
    vector<double> totalParalyticCases;     // vector for num paralytic cases
    vector<double> pCasesPerYear;           // vector for paralytic cases per year
    vector<double> histogramCases(50,0);    // vector for counting number of cases per year

    //int seed = 0;
    //mt19937 gen(seed);

    random_device rd;                       // generates a random real number for the seed
    mt19937 gen(rd());                      // random number generator

    const map<string, double> compartments = initialize_compartments();
    
    const int S  = compartments.at("S");
    const int I1 = compartments.at("I1");
    const int R  = compartments.at("R");
    const int P  = compartments.at("P");
    const int Ir = compartments.at("Ir");

    //The Simulation
    for(int i=0;i<numSims;++i){
        //reset all parameters to original values after each run of the simulation
        double S1=S;
        double I11=I1;
        double R1=R;
        double P1=P;
        double Ir1=Ir;
        double tsc=0;
        double time=0;
        int countPIR = 0;
        pCaseDetection.clear();
        pCasesPerYear.clear();
        
        //run the simulation for 1 mill steps
        for(int j=0;j<1000001;++j){
            double totalRate = 0;
            EventType event_type = sample_event(gen, totalRate, S1, I11, R1, P1, Ir1);
            //Pick the event that is to occur based on generated number and rate
            switch (event_type) {
                case FIRST_INFECTION_EVENT: --S1; ++I11;
                    //generate a unif random real num to determine if a paralytic case is detected
                    {
                        double rr = unifdis(gen);
                        if(rr<(PIR*DET_RATE)){
                            countPIR++;
                            if(countPIR > 1){
                                pCaseDetection.push_back(time-tsc);
                                //tsc=time;
                            }
                            const int year = (int) time;
                            if((unsigned) year >= pCasesPerYear.size()){
                                pCasesPerYear.resize(year + 1, 0);
                            }
                            pCasesPerYear[year]++;
                            tsc=time;
                        }
                    }
                    break;
                case REINFECTION_EVENT:                    --P1; ++Ir1; break;
                case RECOVERY_FROM_FIRST_INFECTION_EVENT:  --I11; ++R1; break;
                case RECOVERY_FROM_REINFECTION_EVENT:      --Ir1; ++R1; break;
                case WANING_EVENT:                         --R1;  ++P1; break;
                case BIRTH_EVENT:                          ++S1;        break;
                case DEATH_FROM_RECOVERED_EVENT:           --R1;        break;
                case DEATH_FROM_PARTIAL_SUSCEPTIBLE_EVENT: --P1;        break;
                case DEATH_FROM_FULLY_SUSCEPTIBLE_EVENT:   --S1;        break;
                case DEATH_FROM_REINFECTION_EVENT:         --Ir1;       break;
                case DEATH_FROM_FIRST_INFECTION_EVENT:     --I11;       break;
                default:
                    cerr << "ERROR: Unsupported event type" << endl;
                    break;
            }

            //generate the time at which the event occurs
            exponential_distribution<>rng(totalRate);
            time+=rng(gen);

            output_streams[S_OUT] << S1      << ", ";
            output_streams[I1_OUT] << I11     << ", ";
            output_streams[R_OUT] << R1      << ", ";
            output_streams[P_OUT] << P1      << ", ";
            output_streams[IR_OUT] << Ir1     << ", ";
            output_streams[TIME_OUT] << time    << ", ";

            //stopping condition
            if((I11+Ir1==0) or time>15){
                totalParalyticCases.push_back(countPIR);
                TTE.push_back(time);
                const double fractional_year = time - (int) time;
                if (fractional_year > 0) {
                    pCasesPerYear.resize((int) time + 1, 0);
                }
                for (auto count: pCasesPerYear) histogramCases[count]++;
                if (time != (int) time) {
                    const int last_years_count = pCasesPerYear.back();
                    if (last_years_count != 0) {
                        histogramCases[last_years_count] -= 1.0 - fractional_year;
                    } else {
                        histogramCases[0] += fractional_year;
                    }
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
                    for (double ptally: histogramCases) output_streams[PCASE_TALLY_OUT] << ptally << ","; output_streams[PCASE_TALLY_OUT] << endl;
                }
                
                for(unsigned int i = 0; i < pCaseDetection.size(); i++){
                    output_streams[PCASE_INTERVAL_OUT] << pCaseDetection[i] << " , ";
                }
                output_streams[PCASE_INTERVAL_OUT] << "\n";
                output_streams[S_OUT] << S1 <<",\n";
                output_streams[I1_OUT] << I11 << ", \n ";
                output_streams[R_OUT] << R1 << ", \n ";
                output_streams[P_OUT] << P1 << ", \n ";
                output_streams[IR_OUT] << Ir1 << ", \n ";
                output_streams[TIME_OUT] << time << ", \n";
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
