#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
    assert(argc == 2);
    string circ_ivl_filename = argv[1];
    vector<double> circ_intervals;
    vector<double> intercase_intervals;
    int num_replicates = 0;
    int required_pcases = 0;
    string ed_stat_output_filename = "e_and_d_values-2pcase_filter-alt-cpp2.out";

    string line;
    ifstream circ_fh(circ_ivl_filename);
    double val;
    if (circ_fh) {
        while (getline(circ_fh, line)) {
            vector<double> times;
            stringstream ss(line);
            while (ss >> val) {
                times.push_back(val);
                if (ss.peek() == ',') ss.ignore();
            }
            vector<double> all_intervals;
            for (unsigned int i = 1; i < times.size(); ++i) {
                all_intervals.push_back(times[i] - times[i-1]); 
            }
            vector<double> intervals;
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

    int num_ci = circ_intervals.size();
    int num_ii = intercase_intervals.size();

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
