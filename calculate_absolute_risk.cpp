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
    string output_filename = "absolute_risk_matrix_N_5x6000_mig_0_fast_v_slow.csv";
    
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
    
    vector<vector<double>> absRiskMat(onePcase_ED_stat.size(),vector<double> (2));
    //vector<vector<double>> absRiskMat(noPcase_ED_stat.size(),vector<double> (2));
    
    ofstream fo(output_filename);
    
    for(int i = 0; i < noPcase_ED_stat.size(); i++){
        for(int j = 0; j < onePcase_ED_stat.size(); j++){
            if(i > 1){
                if(onePcase_ED_stat[j][0] <= noPcase_ED_stat[i][0] and onePcase_ED_stat[j][0]> noPcase_ED_stat[(i-1)][0]){
                    absRiskMat[j][0] = onePcase_ED_stat[j][0];
                    absRiskMat[j][1] = onePcase_ED_stat[j][1] - noPcase_ED_stat[i][1];
                }
            }
            else{
                if(onePcase_ED_stat[j][0] <= noPcase_ED_stat[i][0]){
                    absRiskMat[j][0] = onePcase_ED_stat[j][0];
                    absRiskMat[j][1] = onePcase_ED_stat[j][1] - noPcase_ED_stat[i][1];
                }
            }
        }
    }
    
    /*for(int i = 0; i < onePcase_ED_stat.size(); i++){
        for(int j = 0; j < noPcase_ED_stat.size(); j++){
            if(i > 1){
                if(noPcase_ED_stat[j][0] <= onePcase_ED_stat[i][0] and noPcase_ED_stat[j][0]> onePcase_ED_stat[(i-1)][0]){
                    absRiskMat[j][0] = noPcase_ED_stat[j][0];
                    absRiskMat[j][1] = onePcase_ED_stat[i][1] - noPcase_ED_stat[j][1];
                }
            }
            else{
                if(noPcase_ED_stat[j][0] <= onePcase_ED_stat[i][0]){
                    absRiskMat[j][0] = noPcase_ED_stat[j][0];
                    absRiskMat[j][1] = onePcase_ED_stat[i][1] - noPcase_ED_stat[j][1];
                }
            }
        }
    }*/
    
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


