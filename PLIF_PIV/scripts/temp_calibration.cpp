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

const int total_images = 1000;

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

void linearRegression(vector<double> const x, vector<double> const y, vector<double> &param) {
    //Средние значения
    double Sx = 0, Sy = 0, Sxy = 0, Sxx = 0;
    
    for(int i = 0; i < x.size(); i++){
        Sx += x[i];
        Sy += y[i];
        Sxy += x[i] * y[i];
        Sxx += x[i] * x[i];
    }
    Sx /= x.size();
    Sy /= x.size();
    Sxy /= x.size();
    Sxx /= x.size();
    
    param[0] = (Sx * Sy - Sxy) / (Sx * Sx - Sxx);
    param[1] = (Sxy - param[0] * Sxx) / Sx;
    
}

void createFolders() {
    boost::filesystem::create_directories("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_n/");
}

void getFolder(int img_number, std::string &img_name, std::string &output_temp_image_name) {
    std::string type = ".tif";
    std::string input;
    std::string output;
    
    stringstream str_in;
    stringstream str_out;
    
    if (img_number < 10) {
        input = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100000";
        output = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_n/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100000";
    }
    else if (img_number < 100) {
        input = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S00010000";
        output = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_n/S1_Re5000_43A_T800_2D_25C_p1_C001H001S00010000";
    }
    else if (img_number < 1000) {
        input = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001000";
        output = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_n/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001000";
    }
    else if (img_number < 10000) {
        input = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100";
        output = "/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_n/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100";
    }
    
    str_in << input << img_number << type;
    str_out << output << img_number << type;
    img_name = str_in.str();
    output_temp_image_name = str_out.str();
    str_in.str("");
    str_out.str("");
    
}

long double averagePhoto(Mat const &image) {
    long double average = 0;
    
    for (int i = 0; i < image.rows; i++) {
        for (int j = 0; j < image.cols; j++) {
            average += image.at<ushort>(i,j);
        }
    }
    
    return average;
}

vector<vector<int>> findBlackDots(Mat const image, int temp_treshold) {
    vector<vector<int>> black_dots_coord;
    int k = 0;
    
    for (int i = 0; i < image.rows; i++) {
        for (int j = 0; j < image.cols; j++) {
            if ((image.at<ushort>(i, j) < temp_treshold) || (image.at<ushort>(i, j) > 10 * temp_treshold)) {
                black_dots_coord.push_back(vector<int>(2, 0));
                black_dots_coord[k][0] = i;
                black_dots_coord[k][1] = j;
                k++;
            }
        }
    }
    return black_dots_coord;
}

void dotsImage(Mat const image, cv::String const path, int temp_treshold) {
    Mat tempPhoto = Mat::zeros(488, 1024, CV_8U);
    
    for (int i = 0; i < image.rows; i++) {
        for (int j = 0; j < image.cols; j++) {
            if (image.at<float>(i, j) < temp_treshold)
                tempPhoto.at<uchar>(i,j) = 255;
            else
                tempPhoto.at<uchar>(i,j) = 0;
        }
    }
    
    imwrite(path, tempPhoto);
}

void tempCalibration() {
    ifstream input("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/project/input/temperatures_more.txt");
    cv::String path("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/project/input/average_calibration_croped_2D/*.tif"); //выбор только tif
    vector<cv::String> fn;
    vector<Mat> tempPhoto;
    glob(path,fn,true);
    
    for (size_t k = 0; k < fn.size(); k++) {
        Mat im = cv::imread(fn[k], CV_LOAD_IMAGE_ANYDEPTH);
        
        if (im.empty()) continue; //Пропуск пустых изображений
        tempPhoto.push_back(im);
        
        if (tempPhoto[k].rows == 0) {
            cout << "0 MATRIX!! AVG" << endl;
            return 0;
        }
    }
    
    ////////////////////////////////////
    Mat result_temp_image(878, 1024, CV_32F);
    Mat result_temp_image_int(878, 1024, CV_16U);
    Mat tempPhoto4 = Mat::zeros(786, 1024, CV_32F);
    Mat tempPhoto5 = Mat::zeros(786, 1024, CV_32F);
    Mat tempPhoto5_ = Mat::zeros(92, 1024, CV_32F);
    /////////////////////////////////////
    
    vector<double> temp; //vector
    string str;
    while (getline(input, str)) {
        double _tempr;
        istringstream iss(str);
        iss >> _tempr;
        if (_tempr == 0)
            continue;
        temp.push_back(_tempr);
    }
    input.close();
    
    if (tempPhoto.size() != temp.size())
        cout << "ERROR! CHECK COUNT OF TEMPERATURES AND PHOTOS!" << endl;
    
    vector<double> f_temp = temp;
    vector<double> param(2);
    param[0] = 0;
    param[1] = 0;
    double temp_a = 0;
    double average_ = 0;
    double av_temp = 0;
    int img_number = 1;
    std::string img_name;
    std::string output_temp_image_name;
    std::string output_temp_txt_name;
    
    Rect ROI(0, 92, 1024, 786);
    
    for(int i = 0; i < tempPhoto.size(); i++)
        tempPhoto[i] = tempPhoto[i](ROI);
    
    for(int i = 0; i < temp.size(); i++)
        av_temp += temp[i];
    
    av_temp /= temp.size();
    
    for(int i = 0; i < tempPhoto.size(); i++) //tempPhoto.size!!
        temp[i] = averagePhoto(tempPhoto[i]);
    
    for (int i = 0; i < tempPhoto4.rows; i++) {
        for (int j = 0; j < tempPhoto4.cols; j++) {
            for(int k = 0; k < tempPhoto.size(); k++)
                tempPhoto4.at<ushort>(i,j) += tempPhoto[k].at<ushort>(i,j);
            
            tempPhoto4.at<ushort>(i,j) /= tempPhoto.size();
        }
    }
    
    temp_a = averagePhoto(tempPhoto4);
    
    for(int i = 0; i < tempPhoto.size(); i++) //tempPhoto.size!!
        temp[i] /= temp_a;
    
    linearRegression(f_temp, temp, param);
    
    while (img_number <= total_images) {
        getFolder(img_number, img_name, output_temp_image_name);
        Mat calibratedPhoto = imread(img_name, CV_LOAD_IMAGE_ANYDEPTH);
        calibratedPhoto = calibratedPhoto(ROI);
        
        double t = param[0] * av_temp + param[1]; //Должно быть около 1 (нужно для проверки)
        
        for (int i = 0; i < calibratedPhoto.rows; i++) {
            for (int j = 0; j < calibratedPhoto.cols; j++) {
                tempPhoto5.at<float>(i, j) = static_cast<double>(calibratedPhoto.at<ushort>(i, j)) / tempPhoto4.at<ushort>(i, j);
                tempPhoto5.at<float>(i, j) = (tempPhoto5.at<float>(i, j) - param[1]) / param[0];
            }
        }
        
        for (int i = 0; i < tempPhoto5.rows; i++) {
            for (int j = 0; j < tempPhoto5.cols; j++) {
                average_ += tempPhoto5.at<float>(i,j);
            }
        }
        average_ /= 1024 * 786;
        
        for (int i = 0; i < tempPhoto5_.rows; i++) {
            for (int j = 0; j < tempPhoto5_.cols; j++) {
                tempPhoto5_.at<float>(i,j) = average_;
            }
        }
        
        vconcat(tempPhoto5_, tempPhoto5, result_temp_image);
        
        double atem = 0;
        unsigned short _atem;
        for (int i = 0; i < result_temp_image_int.rows; i++) {
            for (int j = 0; j < result_temp_image_int.cols; j++) {
                atem = 100 * result_temp_image.at<float>(i,j);
                _atem = static_cast<unsigned short>(atem);
                result_temp_image_int.at<ushort>(i,j) = _atem;
            }
        }
        
        imwrite(output_temp_image_name, result_temp_image_int);
        
        if (img_number % 10)
            cout << img_number * 100 / total_images << " %..." << endl;
        
        img_number++;
        average_ = 0;
    }
}

void removeParticles() {
    cv::String path("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature/*.tif");
    cout << "              removePrticles function starts              " << endl;
    vector<cv::String> fn;
    vector<Mat> tempPhoto;
    glob(path, fn, true);
    
    for (size_t k = 0; k < fn.size(); k++)
    {
        Mat im = cv::imread(fn[k], CV_LOAD_IMAGE_ANYDEPTH);
        if (im.empty()) continue;
        tempPhoto.push_back(im);
    }
    
    cout << "              All images are loaded!      " << endl << "              vector size: " << tempPhoto.size() << endl;
    
    for (int i = 0; i < fn.size(); i++) {
        std::string type = ".tif";
        std::string output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_without_dots/S0P41_Re5000_43A_T800_2D_25C_p1_C001H001S000100";
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
        
        vector<vector<int>> black_dots_coord = findBlackDots(tempPhoto[i], 2500);
        
        for (int j = 0; j < black_dots_coord.size(); j++) {
            if (i < fn.size() - 6) {
                tempPhoto[i].at<ushort>(black_dots_coord[j][0], black_dots_coord[j][1]) = tempPhoto[i + 6].at<ushort>(black_dots_coord[j][0], black_dots_coord[j][1]);
            }
            else{
                tempPhoto[i].at<ushort>(black_dots_coord[j][0], black_dots_coord[j][1]) = tempPhoto[i - 6].at<ushort>(black_dots_coord[j][0], black_dots_coord[j][1]);
            }
        }
        
        vector<vector<int>> black_dots_coord2 = findBlackDots(tempPhoto[i], 2500);
        for (int j = 0; j < black_dots_coord2.size(); j++) {
            if (i < fn.size() - 20)
                tempPhoto[i].at<ushort>(black_dots_coord2[j][0], black_dots_coord2[j][1]) = tempPhoto[i + 20].at<ushort>(black_dots_coord2[j][0], black_dots_coord2[j][1]);
            else
                tempPhoto[i].at<ushort>(black_dots_coord2[j][0], black_dots_coord2[j][1]) = tempPhoto[i - 20].at<ushort>(black_dots_coord2[j][0], black_dots_coord2[j][1]);
        }
        
        if (i % 10)
            cout << i * 100 / fn.size() << " %..." << endl;
        imwrite(img_name_1, tempPhoto[i]);
    }
    
}

void fftCalibration() {
    cv::String path("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_without_dots/*.tif");
    cout << "              FFT calibration starts              " << endl;
    vector<cv::String> fn;
    vector<Mat> tempPhoto;
    glob(path, fn, true);
    
    for (size_t k = 0; k < fn.size(); k++) {
        Mat im = cv::imread(fn[k], CV_LOAD_IMAGE_ANYDEPTH);
        if (im.empty()) continue;
        tempPhoto.push_back(im);
    }
    cout << "              All images are loaded!      " << endl << "              vector size: " << tempPhoto.size() << endl;
    
    bool REAL = 1, IMAG = 0;
    fftw_complex x[fn.size()];
    fftw_complex y[fn.size()];
    fftw_complex z[fn.size()];
    
    for (int k = 0; k < tempPhoto[0].rows; k++) {
        for (int j = 0; j < tempPhoto[0].cols; j++) {
            for (int i = 0; i < fn.size(); i++) {
                x[i][REAL] = tempPhoto[i].at<ushort>(k, j);
                x[i][IMAG] = 0;
            }
            
            fftw_plan plan = fftw_plan_dft_1d(static_cast<int>(fn.size()), x, y, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
            fftw_cleanup();
            
            for (int i = 0; i < fn.size(); i++) {
                
                if (i > fn.size() * 0.15) {
                    y[i][REAL] = 0;
                    y[i][IMAG] = 0;
                }
                
            }
            
            fftw_plan plan2 = fftw_plan_dft_1d(static_cast<int>(fn.size()), y, z, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(plan2);
            fftw_destroy_plan(plan2);
            fftw_cleanup();
            
            for (int i = 0; i < fn.size(); i++) {
                tempPhoto[i].at<ushort>(k, j) = static_cast<ushort>(z[i][REAL] / fn.size());
            }
        }
        if (k % 10)
            cout << k * 100.0 / tempPhoto[0].rows << " %..." << endl;
    }
    
    for (int i = 0; i < fn.size(); i++) {
        
        std::string type = ".tif";
        std::string output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_calibration/S0P41_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_without_dots_fft/S0P41_Re5000_43A_T800_2D_25C_p1_C001H001S000100";
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
    tempCalibration();
    removeParticles();
    fftCalibration();
    
    return 0;
}


/*
 Mat tempPhoto0 = imread("/Users/MacBook/Desktop/input/Calibration_curve/AVG_Re5000_43A_noT_2D_Chan_24p70_Box_25p04_p2_C001H001S0P41001.tif", CV_LOAD_IMAGE_ANYDEPTH);
 Mat tempPhoto1 = imread("/Users/MacBook/Desktop/input/Calibration_curve/AVG_Re5000_43A_noT_2D_Chan_25p75_Box_26p03_p2_C001H001S0P41001.tif", CV_LOAD_IMAGE_ANYDEPTH);
 Mat tempPhoto2 = imread("/Users/MacBook/Desktop/input/Calibration_curve/AVG_Re5000_43A_noT_2D_Chan_28p04_Box_28p26_p2_C001H001S0P41P41001.tif", CV_LOAD_IMAGE_ANYDEPTH);
 Mat tempPhoto3 = imread("/Users/MacBook/Desktop/input/Calibration_curve/AVG_Re5000_43A_noT_2D_Chan_29p95_Box_30p05_p3_C001H001S0P41001.tif", CV_LOAD_IMAGE_ANYDEPTH);
 */


/*
 for (int i = 0; i < tempPhoto4.rows; i++) {
 for (int j = 0; j < tempPhoto4.cols; j++) {
 tempPhoto5.at<float>(i, j) = t * tempPhoto4.at<ushort>(i, j);
 }
 }
 
 imwrite("/Users/MacBook/Desktop/output/CHECK.tif", tempPhoto5);
 
 for (int i = 0; i < tempPhoto0.rows; i++) {
 for (int j = 0; j < tempPhoto0.cols; j++) {
 tempPhoto5.at<float>(i, j) = static_cast<double>(tempPhoto0.at<ushort>(i, j)) / tempPhoto4.at<ushort>(i, j);
 tempPhoto5.at<float>(i, j) = (tempPhoto5.at<float>(i, j) - approx_b) / approx_a;
 }
 }
 
 imwrite("/Users/MacBook/Desktop/output/CHECK0.tif", tempPhoto5);
 
 for (int i = 0; i < tempPhoto1.rows; i++) {
 for (int j = 0; j < tempPhoto1.cols; j++) {
 tempPhoto5.at<float>(i, j) = static_cast<double>(tempPhoto1.at<ushort>(i, j)) / tempPhoto4.at<ushort>(i, j);
 tempPhoto5.at<float>(i, j) = (tempPhoto5.at<float>(i, j) - approx_b) / approx_a;
 }
 }
 
 imwrite("/Users/MacBook/Desktop/output/CHECK1.tif", tempPhoto5);
 
 for (int i = 0; i < tempPhoto2.rows; i++) {
 for (int j = 0; j < tempPhoto2.cols; j++) {
 tempPhoto5.at<float>(i, j) = static_cast<double>(tempPhoto2.at<ushort>(i, j)) / tempPhoto4.at<ushort>(i, j);
 tempPhoto5.at<float>(i, j) = (tempPhoto5.at<float>(i, j) - approx_b) / approx_a;
 }
 }
 
 imwrite("/Users/MacBook/Desktop/output/CHECK2.tif", tempPhoto5);
 
 for (int i = 0; i < tempPhoto3.rows; i++) {
 for (int j = 0; j < tempPhoto3.cols; j++) {
 tempPhoto5.at<float>(i, j) = static_cast<double>(tempPhoto3.at<ushort>(i, j)) / tempPhoto4.at<ushort>(i, j);
 tempPhoto5.at<float>(i, j) = (tempPhoto5.at<float>(i, j) - approx_b) / approx_a;
 }
 }
 
 imwrite("/Users/MacBook/Desktop/output/CHECK3.tif", tempPhoto5);
 
 void checkSamePixels(double* f_temp) {
 
 for (int i = 0; i < 4; i++) {
 if (f_temp[i] < f_temp[i + 1])
 f_temp[i] += (f_temp[i + 1] - f_temp[i]) * 2;
 }
 
 for (int i = 0; i < 3; i++) {
 for (int j = i + 1; j < 4; j++) {
 if ((f_temp[i] == f_temp[j]))
 f_temp[i] += 5;
 else if (f_temp[i] == f_temp[j] + 1)
 f_temp[i] += 4;
 else if (f_temp[i] == f_temp[j] + 2)
 f_temp[i] += 3;
 }
 }
 }
 
 //cout << result_temp_image.at<float>(250, 250) << endl;
 //cout << result_temp_image.at<float>(250, 1000) << endl;
 //cout << black_dots_coord.size() << endl;
 //cout << black_dots_coord2.size() << endl;
 
 
 */
