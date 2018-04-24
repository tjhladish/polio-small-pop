#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include <array>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

string output_dir = "/home/tjhladish/work/polio-small-pop/output/";

//const vector<double> kappas = {0.1, 0.2, 0.3};
//const vector<double> rhos   = {10, 20, 30};

int main(int argv, char** argc){
    //waning parameters
    const double kappa = 0.4179; //kappas[atoi(argv[1])];
    const int rho = 0.2; //rhos[atoi(argv[1])];

    //other parameters
    const int recovery = 13;//gamma
    const int beta =135;
    const double birth = 0.02;
    const double death = 0.02;
    //double PIR =0.005; //type 1 paralysis rate
    //const int alpha = 0;

    //initial population from equilibrium values
    const int Tot = 10000;
    const int S =2;
    const int I1=15;
    const int R =2825;
    const int P =649;
    const int Ir=6509;

    //counts

    double AVG_numBirths=0;
    double AVG_numDeaths=0;
    double AVG_numI1Inf=0;
    double AVG_numIrInf=0;
    double AVG_numRec=0;
    double AVG_numWane=0;

    //Number of Simulations to run:
    const int numSims=1000;


    //array holds count of number of first infected (I1) individuals
    array<double,10000> AVG_I1vec = {0};
    AVG_I1vec[0]=I1;

    //array holds count of number of reinfected (Ir) individuals
    array<double,10000> AVG_Irvec = {0};
    AVG_Irvec[0]=Ir;

    //array holds count of number of susceptible (S) individuals
    array<double,10000> AVG_Svec = {0};
    AVG_Svec[0]=S;

    //array holds count of number of recoverd (R) individuals
    array<double,10000> AVG_Rvec = {0};
    AVG_Rvec[0]=R;

    //array holds count of number of partially susceptible (P) indiviudals
    array<double,10000> AVG_Pvec={0};
    AVG_Pvec[0]=P;


    //generate a random real number
    random_device rd;//this is the seed
    //int rd =0;
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
        //double tsc=0;
        double time=0;
        //double z=0;
        //int y=1;
        //int k=1;
        int numBirths=0;
        int numDeaths=0;
        int numI1Inf=0;
        int numIrInf=0;
        int numRec=0;
        int numWane=0;
        AVG_I1vec[0] = I1;
        AVG_Irvec[0] = Ir;
        AVG_Svec[0] = S;
        AVG_Rvec[0] = R;
        AVG_Pvec[0] = P;


        //run the simulation for 1 mill steps
        for(int j=0;j<1000001;++j){

            //update the transition rates
            double birthRate;
            if((S1+I11+R1+P1+Ir1<Tot)){
                birthRate=birth*(S1+I11+R1+P1+Ir1);
            }
            else{
                birthRate=0;
            }
            double infect1 = beta*S1/Tot*(I11+ kappa*Ir1);
            double recover1 = recovery*I11;
            double wane = rho*R1;
            double infectr = kappa*beta*P1/Tot*(I11+ kappa*Ir1);
            double recover2 = (recovery/kappa)*Ir1;
            double deathS=death*S1;
            double deathI1=death*I11;
            double deathR=death*R1;
            double deathP=death*P1;
            double deathIr=death*Ir1;

            double totalRate = birthRate+infect1+recover1+wane+infectr+recover2+deathS+deathI1+deathR+deathP+deathIr;


            //generate unifrn
            double ran=unifdis(gen);

            //Pick the event that is to occur based on generated number and rate
            if(ran< infect1/totalRate){
                S1=S1-1;
                I11=I11+1;
                numI1Inf++;
            }
            else if(ran<((infectr+infect1)/totalRate)){
                P1=P1-1;
                Ir1=Ir1+1;
                numIrInf++;
            }
            else if(ran<((recover1+infectr+infect1)/totalRate)){
                I11=I11-1;
                R1=R1+1;
                numRec++;
            }
            else if(ran<((recover2+recover1+infectr+infect1)/totalRate)){
                Ir1=Ir1-1;
                R1=R1+1;
                numRec++;
            }
            else if(ran<((wane+recover2+recover1+infectr+infect1)/totalRate)){
                R1=R1-1;
                P1=P1+1;
                numWane++;
            }
            else if(ran<((birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                S1=S1+1;
                numBirths++;
            }
            else if(ran<((deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                R1=R1-1;
                numDeaths++;
            }
            else if(ran<((deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                P1=P1-1;
                numDeaths++;
            }
            else if(ran<((deathS+deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                S1=S1-1;
                numDeaths++;
            }
            else if(ran<((deathIr+deathS+deathP+deathR+birthRate+wane+recover2+recover1+infectr+infect1)/totalRate)){
                Ir1=Ir1-1;
                numDeaths++;
            }
            else{
                I11=I11-1;
                numDeaths++;
            }

            //once the event is chosen, generate the time at which the event occurs
            exponential_distribution<>rng((totalRate));

            double t1 = rng(gen);
            time+=t1;

            //stopping condition
            if((I11+Ir1==0) or time>.1){
                AVG_numBirths = (AVG_numBirths*i+numBirths)/(double)(i+1);
                AVG_numDeaths = (AVG_numDeaths*i+numDeaths)/(double)(i+1);
                AVG_numI1Inf = (AVG_numI1Inf*i+numI1Inf)/(double)(i+1);
                AVG_numIrInf = (AVG_numIrInf*i+numIrInf)/(double)(i+1);
                AVG_numRec = (AVG_numRec*i+numRec)/(double)(i+1);
                AVG_numWane= (AVG_numWane*i+numWane)/(double)(i+1);
                break;

            }
        }

    }

    cout<<"num births "<<AVG_numBirths<<"\n";
    cout<<"num deaths "<<AVG_numDeaths<<"\n";
    cout<<"num I1 inf "<<AVG_numI1Inf<<"\n";
    cout<<"num Ir inf "<<AVG_numIrInf<<"\n";
    cout<<"num rec "<<AVG_numRec<<"\n";
    cout<<"num wane "<<AVG_numWane<<"\n";
    ofstream myfile7;
    myfile7.open (output_dir + "N_10000,beta_135,fast_I1vec_time_15_birth_.02_valid.csv");
    for(unsigned int i=0;i<AVG_I1vec.size();++i){
        myfile7<<AVG_I1vec[i]<<"\n";
    }
    myfile7.close();
    ofstream myfile8;
    myfile8.open (output_dir + "N_10000,beta_135,fast_Irvec_time_15_birth_.02_valid.csv");
    for(unsigned int i=0;i<AVG_Irvec.size();++i){
        myfile8<<AVG_Irvec[i]<<"\n";
    }
    myfile8.close();

    ofstream myfile4;
    myfile4.open (output_dir + "N_10000,beta_135,fast_Svec_test_time_15_birth_.02_valid.csv");
    for(unsigned int i=0;i<AVG_Svec.size();++i){
        myfile4<<AVG_Svec[i]<<"\n";
    }
    myfile4.close();
    ofstream myfile5;
    myfile5.open (output_dir + "N_10000,beta_135,fast_Rvec_test_time_15_birth_.02_valid.csv");
    for(unsigned int i=0;i<AVG_Rvec.size();++i){
        myfile5<<AVG_Rvec[i]<<"\n";
    }
    myfile5.close();

    ofstream myfile6;
    myfile6.open (output_dir + "N_10000,beta_135,fast_Pvec_test_time_15_birth_.02_valid.csv");
    for(unsigned int i=0;i<AVG_Pvec.size();++i){
        myfile6<<AVG_Pvec[i]<<"\n";
    }
    myfile6.close();

    return 0;
}
