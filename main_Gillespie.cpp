#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include <array>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <gsl/gsl_linalg.h>

using namespace std;

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
string output_dir ="/Users/Celeste/Desktop/C++PolioSimResults/";

const vector<double> kappas = {0.4179, 0.6383, 0.8434};//fast, intermed, slow
const vector<double> rhos   = {0.2, 0.04, 0.02};//fast, intermed, slow

int main(){
    ofstream myfile;
    ofstream myfile1;
    ofstream myfile2;
    //waning parameters
    const double kappa = 0.4179; //kappas[atoi(argv[1])];
    const double rho = 0.2; //rhos[atoi(argv[1])];

    //other parameters
    const int recovery = 13;//gamma
    const int beta = 135;
    const double birth = 0.02;
    const double death = 0.02;
    double PIR =0.005; //type 1 paralysis rate
    double detRate = 1;
    //const int alpha = 0;

    //initial population from equilibrium values
    const int Tot = 10000;
    const int S =579;
    const int I1=14;
    const int R =4111;
    const int P =5273;
    const int Ir=23;

    myfile.open(output_dir + "time_between_pcases_N_"+to_string(Tot)+",beta_"+to_string(beta)+",detect_rate_"+to_string(detRate)+"rho_"+to_string(rho)+".csv");
    //Number of Simulations to run:
    const int numSims=1000;
    
    //time to extinction vector
    vector<double> TTE;
    
    //time between paralytic cases vector
    vector<double> pCaseDetection;
    
    //vector for num paralytic cases
    vector<double> totalParalyticCases;

    //generate a random real number
    random_device rd;//this is the seed
    //int rd =0;
    //mt19937 gen(1);
    mt19937 gen(rd());//this generates the rn with the above seed
    uniform_real_distribution<> unifdis(0.0, 1.0);

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
        


        //run the simulation for 1 mill steps
        for(int j=0;j<1000001;++j){

            //update the transition rates
            double birthRate;
            if((S1+I11+R1+P1+Ir1<Tot)){
                birthRate = birth*(S1+I11+R1+P1+Ir1);
            }
            else{
                birthRate = 0;
            }
            double infect1 = beta*S1/Tot*(I11+ kappa*Ir1);
            double recover1 = recovery*I11;
            double wane = rho*R1;
            double infectr = kappa*beta*P1/Tot*(I11+ kappa*Ir1);
            double recover2 = (recovery/kappa)*Ir1;
            double deathS = death*S1;
            double deathI1 = death*I11;
            double deathR = death*R1;
            double deathP = death*P1;
            double deathIr = death*Ir1;

            double totalRate = birthRate+infect1+recover1+wane+infectr+recover2+deathS+deathI1+deathR+deathP+deathIr;


            //generate unifrn
            double ran=unifdis(gen);

            //Pick the event that is to occur based on generated number and rate
            if(ran< infect1/totalRate){
                S1=S1-1;
                I11=I11+1;
                //generate a unif random real num to determine if a paralytic case is detected
                double rr = unifdis(gen);
                if(rr<(PIR*detRate)){
                    countPIR++;
                    if(countPIR > 1){
                        pCaseDetection.push_back(time-tsc);
                        tsc=time;
                    }
                }
            }
            else if(ran<((infectr+infect1)/totalRate)){
                P1=P1-1;
                Ir1=Ir1+1;
            }
            else if(ran<((recover1+infectr+infect1)/totalRate)){
                I11=I11-1;
                R1=R1+1;
            }
            else if(ran<((recover2+recover1+infectr+infect1)/totalRate)){
                Ir1=Ir1-1;
                R1=R1+1;
            }
            else if(ran<((wane+recover2+recover1+infectr+infect1)/totalRate)){
                R1=R1-1;
                P1=P1+1;
            }
            else if(ran<((birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                S1=S1+1;
            }
            else if(ran<((deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                R1=R1-1;
            }
            else if(ran<((deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                P1=P1-1;
            }
            else if(ran<((deathS+deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                S1=S1-1;
            }
            else if(ran<((deathIr+deathS+deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                Ir1=Ir1-1;
            }
            else{
                I11=I11-1;
            }

            //once the event is chosen, generate the time at which the event occurs
            exponential_distribution<>rng((totalRate));

            double t1 = rng(gen);
            time+=t1;

            //stopping condition
            if((I11+Ir1==0) or time>15){
                totalParalyticCases.push_back(countPIR);
                TTE.push_back(time);
                for(unsigned int i = 0; i < pCaseDetection.size(); i++){
                    myfile << pCaseDetection[i] << " , ";
                }
                myfile << "\n";
                break;
            }
        }

    }
    myfile1.open(output_dir + "num_p_cases_N_"+to_string(Tot)+",beta_"+to_string(beta)+",detect_rate_"+to_string(detRate)+"rho_"+to_string(rho)+".csv");
    for(unsigned int i = 0; i < totalParalyticCases.size(); i++){
        myfile1<<totalParalyticCases[i]<<"\n";
    }
    myfile1.close();
    myfile2.open(output_dir + "TTE_N_"+to_string(Tot)+",beta_"+to_string(beta)+",detect_rate_"+to_string(detRate)+"rho_"+to_string(rho)+".csv");
    for (unsigned int i = 0; i < TTE.size(); i++) {
        myfile2<<TTE[i]<<"\n";
    }
    myfile2.close();
    myfile.close();
    return 0;
}
