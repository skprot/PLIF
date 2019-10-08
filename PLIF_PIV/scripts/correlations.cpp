#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#include <boost/filesystem.hpp>


using namespace std;

const int total_images = 1000;

void createFolders() {
    boost::filesystem::create_directories("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/correlations/");
}

string getNameAndFolder(int img_number, short input) {
    
    char name[350];
    const char* folder1;
    const char* folder2;
    const char* folder3;
    const char* folder4;
    
    if (input == 1) {
        folder1 = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_overlaped_noise_txt/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100000%d.txt";
        folder2 = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_overlaped_noise_txt/S1_Re5000_43A_T800_2D_25C_p1_C001H001S00010000%d.txt";
        folder3 = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_overlaped_noise_txt/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001000%d.txt";
        folder4 = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_overlaped_noise_txt/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100%d.txt";
    }
    else if (input == 0) {
        folder1 = "/media/itpuser/IMP_JET_PROC_STOR1/Processing_ImpJet_Heat_PIV/Exp_copys/!Export_Inst_All_v12_18/S1_Re5000_43A_T800_25C_2D/S1_Re5000_2D_000%d.txt";
        folder2 = "/media/itpuser/IMP_JET_PROC_STOR1/Processing_ImpJet_Heat_PIV/Exp_copys/!Export_Inst_All_v12_18/S1_Re5000_43A_T800_25C_2D/S1_Re5000_2D_00%d.txt";
        folder3 = "/media/itpuser/IMP_JET_PROC_STOR1/Processing_ImpJet_Heat_PIV/Exp_copys/!Export_Inst_All_v12_18/S1_Re5000_43A_T800_25C_2D/S1_Re5000_2D_0%d.txt";
        folder4 = "/media/itpuser/IMP_JET_PROC_STOR1/Processing_ImpJet_Heat_PIV/Exp_copys/!Export_Inst_All_v12_18/S1_Re5000_43A_T800_25C_2D/S1_Re5000_2D_%d.txt";
    }
    else if (input == 3) {
        folder1 = "/Users/MacBook/Desktop/output/correlations/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100000%d.txt";
        folder2 = "/Users/MacBook/Desktop/output/correlations/S1_Re5000_43A_T800_2D_25C_p1_C001H001S00010000%d.txt";
        folder3 = "/Users/MacBook/Desktop/output/correlations/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001000%d.txt";
        folder4 = "/Users/MacBook/Desktop/output/correlations/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100%d.txt";
    }
    else {
        cout << "FOLDER ERROR!";
        return 0;
    }
    
    if (img_number < 10) {
        sprintf(name, folder1, img_number);
    }
    else if (img_number < 100) {
        sprintf(name, folder2, img_number);
    }
    else  if (img_number < 1000) {
        sprintf(name, folder3, img_number);
    }
    else {
        sprintf(name, folder4, img_number);
    }
    
    return name;
}

void findCorrelations (vector<vector<double>> plif_temp, vector<vector<double>> plif_x, vector<vector<double>> plif_y, vector<vector<double>> piv_x, vector<vector<double>> piv_y, vector<vector<double>> piv_vx, vector<vector<double>> piv_vy, vector<vector<double>> piv_vz) {
    
    vector<double> correl_1;
    vector<double> correl_2;
    vector<double> correl_3;
    vector<double> correl_4;
    vector<double> correl_5;
    vector<double> correl_6;
    vector<double> correl_7;
    vector<double> correl_8;
    vector<double> correl_9;
    vector<double> correl_10;
    vector<double> correl_11;
    vector<double> vx_T;
    vector<double> vy_T;
    vector<double> vz_T;
    vector<double> avg_vx;
    vector<double> avg_vy;
    vector<double> avg_vz;
    vector<double> avg_vx_vy;
    vector<double> avg_temp;
    vector<double> avg_vx_2;
    vector<double> avg_vy_2;
    vector<double> avg_vz_2;
    vector<double> avg_temp_2;
    
    ofstream output("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/correlations/correl_S1_5450_2D_with_noise.txt");
    output << "[xmm]" << " [ymm]" << " [vx]" << " [vy]" << " [vz]" << " [T]" << " [vx^2]" << " [vy^2]" << " [vz^2]" << " [T^2]" << " [av_T*av_vx]" << " [av_T*av_vy]" << " [av_T*av_vz]" << " [av_(T^2)-av_T^2]" << " [av_T*av_vx-av_vx*T]" << " [av_T*av_vy-av_vy*T]"  << " [av_T*av_vz-av_vz*T]" << " [av_(vx^2)-av_vx^2]" << " [av_(vy^2)-av_vy^2]" << " [av_(vz^2)-av_vz^2]" << " [<vx*vy>-<vx><vy>]" << endl;
    
    for (int it = 0; it < plif_temp[0].size(); it++) {
        double _temp = 0, _temp_2 = 0;
        int iter = 1;
        while(iter <= total_images) {
            _temp = _temp + (plif_temp[iter - 1][it] / 100.0);
            _temp_2 = _temp_2 + (plif_temp[iter - 1][it] * plif_temp[iter - 1][it] / 10000.0);
            iter++;
        }
        
        _temp /= plif_temp.size();
        _temp_2 /= plif_temp.size();
        avg_temp.push_back(_temp);
        avg_temp_2.push_back(_temp_2);
    }
    
    for (int it = 0; it < piv_x[0].size(); it++) { // Все что связано только со скоростями
        double _vx = 0, _vy = 0, _vz = 0, _vx_2 = 0, _vy_2 = 0, _vz_2 = 0, _vx_vy = 0;
        int iter = 1;
        while(iter <= total_images) {
            _vx += piv_vx[iter - 1][it];
            _vy += piv_vy[iter - 1][it];
            _vz += piv_vz[iter - 1][it];
            _vx_2 += piv_vx[iter - 1][it] * piv_vx[iter - 1][it];
            _vy_2 += piv_vy[iter - 1][it] * piv_vy[iter - 1][it];
            _vz_2 += piv_vz[iter - 1][it] * piv_vz[iter - 1][it];
            _vx_vy += piv_vx[iter - 1][it] * piv_vy[iter - 1][it];
            iter++;
        }
        
        _vx /= piv_vx.size();
        _vy /= piv_vy.size();
        _vz /= piv_vz.size();
        _vx_2 /= piv_vx.size();
        _vy_2 /= piv_vx.size();
        _vz_2 /= piv_vx.size();
        _vx_vy /= piv_vx.size();
        
        avg_vx.push_back(_vx);
        avg_vy.push_back(_vy);
        avg_vz.push_back(_vz);
        avg_vx_2.push_back(_vx_2);
        avg_vy_2.push_back(_vy_2);
        avg_vz_2.push_back(_vz_2);
        avg_vx_vy.push_back(_vx_vy);
    }
    
    for (int it = 0; it < plif_temp[0].size(); it++) { // все что связано с температурами и скоростями
        int i = 0;
        while (((abs(piv_x[0][i] - plif_x[0][it])) > 0.05) || (abs(piv_y[0][i] -  plif_y[0][it])) > 0.05) {
            i++;
        }
        
        cout << piv_x[0][i] << " " << plif_x[0][it] << "   " << piv_y[0][i] << " " << plif_y[0][it] << endl;
        correl_1.push_back(avg_temp[it] * avg_vx[i]);
        correl_2.push_back(avg_temp[it] * avg_vy[i]);
        correl_3.push_back(avg_temp[it] * avg_vz[i]);
        correl_4.push_back(avg_temp_2[it] - avg_temp[it] * avg_temp[it]);
        cout << "4" << endl;
        
        int iter = 1;
        double _vx_T = 0, _vy_T = 0, _vz_T = 0;
        while(iter <= total_images) {
            _vx_T += (plif_temp[iter - 1][it] / 100.0) * piv_vx[iter - 1][i];
            _vy_T += (plif_temp[iter - 1][it] / 100.0) * piv_vy[iter - 1][i];
            _vz_T += (plif_temp[iter - 1][it] / 100.0) * piv_vz[iter - 1][i];
            iter++;
        }
        
        vx_T.push_back(_vx_T / (iter - 1));
        vy_T.push_back(_vy_T / (iter - 1));
        vz_T.push_back(_vz_T / (iter - 1));
        
        correl_5.push_back(vx_T[it] - avg_vx[i] * avg_temp[it]);
        correl_6.push_back(vy_T[it] - avg_vy[i] * avg_temp[it]);
        correl_7.push_back(vz_T[it] - avg_vz[i] * avg_temp[it]);
        correl_8.push_back(avg_vx_2[i] - avg_vx[i] * avg_vx[i]);
        correl_9.push_back(avg_vy_2[i] - avg_vy[i] * avg_vy[i]);
        correl_10.push_back(avg_vz_2[i] - avg_vz[i] * avg_vz[i]);
        correl_11.push_back(avg_vx_vy[i] - avg_vx[i] * avg_vy[i]);
        
        output << plif_x[0][it] << " " << plif_y[0][it] << " " << avg_vx[i] << " " << avg_vy[i] << " " << avg_vz[i] << " "  << avg_temp[it] << " "  << avg_vx_2[it] << " "  << avg_vy_2[it]  << " "  << avg_vz_2[it] <<  " "  << avg_temp_2[it] << " "  << correl_1[it] << " " << correl_2[it] << " " << correl_3[it] << " " << correl_4[it] << " " << correl_5[it] << " " << correl_6[it] << " " << correl_7[it]  << " " << correl_8[it]  << " " << correl_9[it]  << " " << correl_10[it] << " " << correl_11[it] << endl;
    }
}

int main() {
    createFolders();
    ifstream infile, piv_infile;
    
    vector<vector<double>> plif_temp(0, vector<double>());
    vector<vector<double>> plif_x(0, vector<double>());
    vector<vector<double>> plif_y(0, vector<double>());
    vector<vector<double>> piv_x(0, vector<double>());
    vector<vector<double>> piv_y(0, vector<double>());
    vector<vector<double>> piv_vx(0, vector<double>());
    vector<vector<double>> piv_vy(0, vector<double>());
    vector<vector<double>> piv_vz(0, vector<double>());
    
    
    int input_lif = 1, input_piv = 0, output = 3;
    
    int num = 0, num_p = 0;
    int iter = 1;
    
    while(iter <= total_images) {
        infile.open(getNameAndFolder(iter, input_lif), ios::in);
        if(infile.fail())
        {
            cout << "error" << endl;
            return 1;
        }
        
        plif_x.push_back({0});
        plif_y.push_back({0});
        plif_temp.push_back({0});
        
        num = 0;
        
        while(!infile.eof())
        {
            double temp = 0, x = 0, y = 0;
            if (num == 0) {
                string temp_s = "";
                getline(infile, temp_s);
                num++;
                continue;
            }
            infile >> x >> y >> temp;
            
            plif_x[iter - 1].push_back(x);
            plif_y[iter - 1].push_back(y);
            if (temp > 3000)
                temp = 2600;
            plif_temp[iter - 1].push_back(temp);
            
            ++num; // go to the next number
        }
        infile.close();
        
        plif_x[iter - 1].erase(plif_x[iter - 1].begin());
        plif_y[iter - 1].erase(plif_y[iter - 1].begin());
        plif_temp[iter - 1].erase(plif_temp[iter - 1].begin());
        
        if (plif_x.size() % 25)
            cout << plif_x.size() * 100.0 / total_images << " %..." << endl;
        
        piv_infile.open(getNameAndFolder(iter, input_piv), ios::in);
        if(piv_infile.fail())
        {
            cout << "error" << endl;
            return 1;
        }
        
        piv_x.push_back({0});
        piv_y.push_back({0});
        piv_vx.push_back({0});
        piv_vy.push_back({0});
        piv_vz.push_back({0});
        
        num_p = 0;
        
        while(!piv_infile.eof())
        {
            double vx = 0, vy = 0, vz = 0, x_p = 0, y_p = 0, tr_1 = 0, tr_2 = 0;
            
            if (num_p < 10) {
                string xp_s = "";
                getline(piv_infile, xp_s);
                num_p++;
                continue;
            }
            
            piv_infile >> x_p >> y_p >> vx >> vy >> vz >> tr_1 >> tr_2;
            
            piv_x[iter - 1].push_back(x_p);
            piv_y[iter - 1].push_back(y_p);
            piv_vx[iter - 1].push_back(vx);
            piv_vy[iter - 1].push_back(vy); //x
            piv_vz[iter - 1].push_back(vz);
            
            
            ++num_p; // go to the next number
            
        }
        piv_infile.close();
        
        piv_x[iter - 1].erase(piv_x[iter - 1].begin());
        piv_y[iter - 1].erase(piv_y[iter - 1].begin());
        piv_vx[iter - 1].erase(piv_vx[iter - 1].begin());
        piv_vy[iter - 1].erase(piv_vy[iter - 1].begin());
        piv_vz[iter - 1].erase(piv_vz[iter - 1].begin());
        piv_x[iter - 1].pop_back();
        piv_y[iter - 1].pop_back();
        piv_vx[iter - 1].pop_back();
        piv_vy[iter - 1].pop_back();
        piv_vz[iter - 1].pop_back();
        plif_x[iter - 1].pop_back();
        plif_y[iter - 1].pop_back();
        plif_temp[iter - 1].pop_back();
        
        iter++;
        
    }
    
    if (abs(piv_x[0][0] + 31.15911) > 0.1) {
        cout << "      first piv coord wrong!:   " << piv_x[0][0] << endl;
        cout << "      size piv coord:   " << piv_x.size() << endl;
        return 0;
    }
    
    findCorrelations(plif_temp, plif_x, plif_y, piv_x, piv_y, piv_vx, piv_vy, piv_vz);
    return 0; // everything went right.
}
