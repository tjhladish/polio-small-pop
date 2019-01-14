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
    if(argc != 3){
        cout<< "argument input: <0pcase file> <1pcase file>"<<endl;
        exit(-1);
    }
    //assert(argc == 3);
    string noPcase_filename = argv[1];
    string onePcase_filename = argv[2];
    
    vector<vector<double>> noPcase_ED_stat;
    vector<vector<double>> onePcase_ED_stat;
    string output_filename = "absolute_risk_matrix_N_5000.csv";
    
    string line;
    double noPcaseValue;
    double onePcaseValue;
    ifstream noPcase_file(noPcase_filename);
    ifstream onePcase_file(onePcase_filename);
    
    while(getline(noPcase_file,line)){
        vector<double> row;
        istringstream iss(line);
        while(iss >> noPcaseValue){
            row.push_back(noPcaseValue);
        }
        noPcase_ED_stat.push_back(row);
    }
    while(getline(onePcase_file,line)){
        vector<double> row;
        istringstream iss(line);
        while(iss >> onePcaseValue){
            row.push_back(onePcaseValue);
        }
        onePcase_ED_stat.push_back(row);
    }
    cout<<"size of no pcase ed stat "<<noPcase_ED_stat.size()<<endl;
    cout<<"size of one pcase ed stat "<<onePcase_ED_stat.size()<<endl;
    
    noPcase_file.close();
    onePcase_file.close();
    
    ofstream fo(output_filename);
    
    vector<double> onePcase_time_vec(onePcase_ED_stat.size());
    vector<double> noPcase_time_vec(noPcase_ED_stat.size());
    
    for(int i = 0; i < onePcase_ED_stat.size(); i++){
        onePcase_time_vec[i] = onePcase_ED_stat[i][0];
    }
    for(int i = 0; i < noPcase_ED_stat.size(); i++){
        noPcase_time_vec[i] = noPcase_ED_stat[i][0];
    }
    
    vector<double> newTimeVec = onePcase_time_vec;
    newTimeVec.insert(newTimeVec.end(), noPcase_time_vec.begin(), noPcase_time_vec.end());
    sort(newTimeVec.begin(),newTimeVec.end());
    
    vector<vector<double>> absRiskMat(newTimeVec.size(),vector<double> (2));
    
    double onePcase_chosen;
    double noPcase_chosen;
    int onePcase_counter = 0;
    int noPcase_counter = 0;
    for(int time = 0; time < newTimeVec.size(); time++){
        onePcase_chosen = 0;
        noPcase_chosen = 0;
        if(onePcase_ED_stat[onePcase_counter][0] > newTimeVec[time]){
            if(onePcase_counter > 0){
                onePcase_chosen = onePcase_ED_stat[(onePcase_counter - 1)][1];
            }
            else{
                onePcase_chosen = onePcase_ED_stat[onePcase_counter][1];
            }
        }
        else{
            onePcase_chosen = onePcase_ED_stat[onePcase_counter][1];
            if(onePcase_counter < (onePcase_ED_stat.size()-1)){
                onePcase_counter++;
            }
        }
        if(noPcase_ED_stat[noPcase_counter][0] > newTimeVec[time]){
            if(noPcase_counter > 0){
                noPcase_chosen = noPcase_ED_stat[(noPcase_counter - 1)][1];
            }
            else{
                noPcase_chosen = noPcase_ED_stat[noPcase_counter][1];
            }
        }
        else{
            noPcase_chosen = noPcase_ED_stat[noPcase_counter][1];
            if(noPcase_counter < (noPcase_ED_stat.size()-1)){
                noPcase_counter++;
            }
        }
        assert(onePcase_chosen!=0);
        assert(noPcase_chosen!=0);
        absRiskMat[time][0] = newTimeVec[time];
        absRiskMat[time][1] = onePcase_chosen - noPcase_chosen;
    }
    
    for(int row = 0; row < absRiskMat.size(); row++){
        for(int col = 0; col < absRiskMat[row].size(); col++){
            fo << absRiskMat[row][col];
            if(col < (absRiskMat[row].size()-1)){
                fo << ", ";
            }
            else{
                fo << endl;
            }
        }
    }
    fo.close();
    return 0;
}


