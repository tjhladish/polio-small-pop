#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <sstream>
#include <random>

using namespace std;

int main(int argc, char** argv) {
    assert(argc == 2);
    string circ_ivl_filename = argv[1];
    vector<double> circ_intervals;
    vector<double> intercase_intervals;
    int num_replicates = 0;
    int required_pcases = 0;
    string ed_stat_output_filename = "e_and_d_values-det_1_N_2x8000_beta_135_fast_migrate_0_birth_with_death.out";
   // "e_and_d_values-0pcase_filter_det_1_N_5000500050005000_beta_135_fast_metapop_4_vil_movement_migrate_0.5.out";
    //"e_and_d_values-0pcase_filter_det_1_N_30000_beta_135_fast_movement_paper.out";
    //exponential_distribution<double> expDis((1/.5359)); //for 3500
    exponential_distribution<double> expDis((1/0.8022)); //for 10,000 E&D first intercase interval
    //exponential_distribution<double> expDis((1/0.7035)); //for 10,000 true intercase
    //exponential_distribution<double> expDis((1/0.6638)); //for 10,000 true intercase
    //exponential_distribution<double> expDis((1/0.9228));//for 10,000 EE
    //exponential_distribution<double> expDis(1.8863);//for 20,000
    //exponential_distribution<double> expDis(.32664); //for 3500
    random_device rd;                       // generates a random real number for the seed
    mt19937 gen(rd());                      // random number generator

    string line;
    ifstream circ_fh(circ_ivl_filename);
    double val;
    if (circ_fh) {
        while (getline(circ_fh, line)) {
            vector<double> times;
            stringstream ss(line);
            while (ss >> val) {
                times.push_back(val);
                if (ss.peek() == ','){
                    ss.ignore();
                }
            }
            vector<double> all_intervals;
            vector<double> intervals;
            for (unsigned int i = 1; i < times.size(); ++i) {
                /*if(i == 1 and times.size()>2){
                    all_intervals.push_back(expDis(gen));
                }*/
                //else{
                    all_intervals.push_back(times[i] - times[i-1]);
                //}
            }
            if (required_pcases < all_intervals.size()) {
                copy(all_intervals.begin()+required_pcases, all_intervals.end(), back_inserter(intervals));
            }
            circ_intervals.insert(circ_intervals.end(), intervals.begin(), intervals.end());
            if (intervals.size() > 0) {
                intercase_intervals.insert(intercase_intervals.end(), intervals.begin(), intervals.end()-1);
                num_replicates++;
            }
        }
    }
    circ_fh.close();


    sort(circ_intervals.begin(), circ_intervals.end()); 
    sort(intercase_intervals.begin(), intercase_intervals.end());
    

    const int num_ci = circ_intervals.size();
    const int num_ii = intercase_intervals.size();

    int last_icase_idx = 0;
    unsigned int icase_idx = 0;


    ofstream fo(ed_stat_output_filename);

    for (unsigned int i = 0; i < num_ci; ++i) {
        double numerator = num_ci - i;
        double interval_duration = circ_intervals[i]; 
        for (icase_idx = last_icase_idx; icase_idx < intercase_intervals.size(); ++icase_idx) {
            if (intercase_intervals[icase_idx] >= interval_duration) {
                break;
            }
        }
        last_icase_idx = icase_idx;
        double denominator = num_replicates + num_ii - icase_idx;
        fo << interval_duration << " " << numerator/denominator << endl;
    }
    fo.close();
}
