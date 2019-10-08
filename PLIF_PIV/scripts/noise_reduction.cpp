#include <iostream>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include <string>
#include <stdint.h>
#include <stdio.h>
#include <sstream>

using namespace cv;
using namespace std;

///////////////////////////////////////////
                                         //
const int total_images = 5459;           //
                                         //
///////////////////////////////////////////

void createFolders() {
    boost::filesystem::create_directories("/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_noise_reduct/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/"); // write here folder for creation
}

void getFolder(int img_number, std::string &img_name, std::string &output_temp_image_name) {
    std::string type = ".tif";
    std::string input;
    std::string output;
    
    stringstream str_in;
    stringstream str_out;
    
    if (img_number < 10) {
        input = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100000";
        output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_noise_reduct/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100000";
    }
    else if (img_number < 100) {
        input = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S00010000";
        output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_noise_reduct/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S00010000";
    }
    else if (img_number < 1000) {
        input = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001000";
        output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_noise_reduct/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001000";
    }
    else if (img_number < 10000) {
        input = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/Series_06_09_17/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100";
        output = "/media/itpuser/IMP_JET_PROC_STOR/Stereo/sayanLIF/Series_06_09_17_noise_reduct/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/S1_Re5000_43A_T800_2D_25C_p1_C001H001S000100";
    }
    
    str_in << input << img_number << type;
    str_out << output << img_number << type;
    img_name = str_in.str();
    output_temp_image_name = str_out.str();
    str_in.str("");
    str_out.str("");
}

void inverseFurierTransformation(Mat& complexI, std::string const &output_temp_image_name) {
    Mat inverseTransform;
    idft(complexI, inverseTransform, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);
    inverseTransform.convertTo(inverseTransform, CV_16U);
    imwrite(output_temp_image_name, inverseTransform);
}

void furierTransformation(Mat& imagePhoto, Mat& mask, std::string const &output_temp_image_name) {
    Mat planes[] = {Mat_<float>(imagePhoto), Mat::zeros(imagePhoto.size(), CV_32F)};
    Mat complexI;
    Mat inverseTransform;
    
    merge(planes, 2, complexI);
    
    dft(complexI, complexI);

    split(complexI, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
 
    for (int i = 0; i < complexI.rows; i++) {
        for (int j = 0; j < complexI.cols; j++) {
            planes[0].at<float>(i,j) *= mask.at<uchar>(i,j);
            planes[1].at<float>(i,j) *= mask.at<uchar>(i,j);
        }
    }
    
    merge(planes, 2, complexI);
    
    idft(complexI, inverseTransform, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);
    inverseTransform.convertTo(inverseTransform, CV_16U);
    imwrite(output_temp_image_name, inverseTransform);
}

int main() {
    Mat image_mask = imread("/home/samba/allaccess/progs/mask_2d.png", CV_LOAD_IMAGE_ANYDEPTH);
    createFolders();
    
    image_mask.convertTo(image_mask, CV_8U);
    
    for (int i = 0; i < image_mask.rows; i++) {
        for (int j = 0; j < image_mask.cols; j++) {
            if (image_mask.at<uchar>(i,j) > 15) {
                image_mask.at<uchar>(i,j) = 1;
            }
            else
                image_mask.at<uchar>(i,j) = 0;
        }
    }
    
    int img_number = 1;
    std::string input, output;
    
    while (img_number <= total_images) {
        getFolder(img_number, input, output);
        Mat image_noise = imread(input, CV_LOAD_IMAGE_ANYDEPTH);
        if (img_number % 30) {
            cout <<  (100.0 / total_images) * img_number  << " %..." << endl;
        }
        
        Rect ROI(0, 0, 1024, 878);
        image_noise = image_noise(ROI);
    
        furierTransformation(image_noise, image_mask, output);
        img_number++;
    }
    return 0;
}

