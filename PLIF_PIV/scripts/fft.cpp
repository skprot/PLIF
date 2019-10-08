#include <iostream>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include <string>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <fftw3.h>

using namespace cv;
using namespace std;

template <class TYPE>
void getPixelValue(Mat& imagePhoto, string name) {
    TYPE pixels[imagePhoto.rows][imagePhoto.cols];
    
    ofstream fout(name);
    
    for (int i = 0; i < imagePhoto.rows; i++) {
        for (int j = 0; j < imagePhoto.cols; j++) {
            pixels[i][j] = imagePhoto.at<TYPE>(i, j);
            fout << pixels[i][j] <<" ";
        }
        fout << "\n";
    }
}


void createFolders() {
    boost::filesystem::create_directories("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_1D_25C_p2_C001H001S0001/temperature_without_dots_fft/");
}

void fftCalibration() {
    cv::String path("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_1D_25C_p2_C001H001S0001/temperature_without_dots/*.tif");

    vector<cv::String> fn;
    vector<Mat> tempPhoto;
    glob(path, fn, true);
    ofstream fout("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/fft.txt");
    ofstream fout1("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/fft1.txt");
    ofstream fout2("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/fft2.txt");

    for (size_t k = 0; k < fn.size(); k++) {
        Mat im = cv::imread(fn[k], CV_LOAD_IMAGE_ANYDEPTH);
        if (im.empty()) continue;
        tempPhoto.push_back(im);
    }
    
    bool REAL = 1, IMAG = 0;
    fftw_complex x[fn.size()];
    fftw_complex y[fn.size()];
    fftw_complex z[fn.size()];

    for (int k = 0; k < tempPhoto[0].rows; k++) {
        cout << "k = " << k << " temp:  " << tempPhoto[10].at<ushort>(k, 250) << endl;
        for (int j = 0; j < tempPhoto[0].cols; j++) {
            for (int i = 0; i < fn.size(); i++) {
                x[i][REAL] = tempPhoto[i].at<ushort>(k, j);
                x[i][IMAG] = 0;
            }

            fftw_plan plan = fftw_plan_dft_1d(static_cast<int>(fn.size()), x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
            fftw_cleanup();
            
            if (j == 250) {
            for (int i = 0; i < fn.size(); i++) {
                fout << y[i][REAL] << endl;
                if (i > fn.size()) {
                    y[i][REAL] = 0;
                    y[i][IMAG] = 0;
                }
                //fout << y[i][REAL] << endl;
            }
            
            for (int i = 0; i < fn.size(); i++) {
                fout1 << y[i][IMAG] << endl;
                if (i > fn.size()) {
                    y[i][REAL] = 0;
                    y[i][IMAG] = 0;
                }
                //fout << y[i][REAL] << endl;
            }
            
            for (int i = 0; i < fn.size(); i++) {
                fout2 << y[i][IMAG]*y[i][IMAG] + y[i][REAL]*y[i][REAL] << endl;
                if (i > fn.size()) {
                    y[i][REAL] = 0;
                    y[i][IMAG] = 0;
                }
                //fout << y[i][REAL] << endl;
            }
            
            fftw_plan plan2 = fftw_plan_dft_1d(static_cast<int>(fn.size()), y, z, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(plan2);
            fftw_destroy_plan(plan2);
            fftw_cleanup();
            
            for (int i = 0; i < fn.size(); i++) {
                tempPhoto[i].at<ushort>(k, j) = static_cast<ushort>(z[i][REAL] / fn.size());
           }
        }
        cout << "k = " << k << " temp:  " << tempPhoto[10].at<ushort>(k, 250) << endl;
    }
    
    for (int i = 0; i < fn.size(); i++) {
        
        std::string type = ".tif";
        std::string output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_1D_25C_p2_C001H001S0001/temperature_without_dots_fft/S0P41_Re5000_43A_T800_1D_25C_p2_C001H001S000100";
        std::string img_name_1;
        
        stringstream str_out;
        
        
        if (i < 9) {
            str_out << output << "000" << i + 1 << type;
        }
        else if (i < 99) {
            str_out << output << "00" << i + 1 << type;
        }
        else if (i < 999) {
            str_out << output << "0" << i + 1 << type;
            
        }
        else {
            str_out << output << i + 1 << type;
        }
        
        img_name_1 = str_out.str();
        str_out.str("");
       
        imwrite(img_name_1, tempPhoto[i]);
    }
}

int main() {
    createFolders();
    fftCalibration();
    
    return 0;
}
