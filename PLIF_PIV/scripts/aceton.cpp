#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <numeric>

using namespace std;
using namespace cv;

const double acenton_conc = 0.0002267; //mol per cm^3


string getNameAndFolder(int img_number, short input) {
    
    char name[124];
    const char* folder1;
    const char* folder2;
    const char* folder3;
    
    if (input == 1) {
        folder1 = "/Users/MacBook/Desktop/input/S1_Re_5000_set2/B0000%d.tif";
        folder2 = "/Users/MacBook/Desktop/input/S1_Re_5000_set2/B000%d.tif";
        folder3 = "/Users/MacBook/Desktop/input/S1_Re_5000_set2/B00%d.tif";
    }
    else if (input == 0) {
        folder1 = "/Users/MacBook/Desktop/output/S1_Re_5000_set2/B0000%d.tiff";
        folder2 = "/Users/MacBook/Desktop/output/S1_Re_5000_set2/B000%d.tiff";
        folder3 = "/Users/MacBook/Desktop/output/S1_Re_5000_set2/B00%d.tiff";
    }
    else if (input == 3) {
        folder1 = "/Users/MacBook/Desktop/output/averaged_with_overlap_images_txt_changed/S1_Re5000_43A_T800_1D_25C_p1_C001H001S000100000%d.txt";
        folder2 = "/Users/MacBook/Desktop/output/averaged_with_overlap_images_txt_changed/S1_Re5000_43A_T800_1D_25C_p1_C001H001S00010000%d.txt";
        folder3 = "/Users/MacBook/Desktop/output/averaged_with_overlap_images_txt_changed/S1_Re5000_43A_T800_1D_25C_p1_C001H001S0001000%d.txt";
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
    else {
        sprintf(name, folder3, img_number);
    }
    
    return name;
}

int main() {
    int img_number = 1;
    bool input = 1, output = 0;
    
    string img_name;
    string output_image_name;
    Mat background = imread("/Users/MacBook/Desktop/input/avg_background/AVG_avg_background.tif", CV_LOAD_IMAGE_ANYDEPTH);
    Mat sheet = imread("/Users/MacBook/Desktop/input/avg_sheet/AVG_sheet_line_1.tif", CV_LOAD_IMAGE_ANYDEPTH);
    Mat aceton = imread("/Users/MacBook/Desktop/input/avg_aceton/AVG_aceton.tif", CV_LOAD_IMAGE_ANYDEPTH);
    Mat aceton_min = imread("/Users/MacBook/Desktop/input/avg_aceton/avg_aceton_min.tif", CV_LOAD_IMAGE_ANYDEPTH);
    Mat new_sheet = Mat::zeros(background.rows, background.cols, CV_16U);
    Mat _new_sheet = Mat::zeros(background.rows, background.cols, CV_32F);


    int j = 0, i = 0;
    
    double av = 0, aceton_value_max = 0, aceton_value_min = 0, k = 0;
    
    for (i = 0; i < sheet.rows; i++) {
        av = 0;
        for (j = 0; j < sheet.cols; j++) {
            av += (sheet.at<ushort>(i,j) - 1300);
        }
        av /= sheet.cols;
        for (j = 0; j < new_sheet.cols; j++) {
            new_sheet.at<ushort>(i,j) = av;
        }
    }
    imwrite("/Users/MacBook/Desktop/output/new_sheet.tif", new_sheet);
    
    for (i = 0; i < aceton.rows; i++) {
        for (j = 0; j < aceton.cols; j++) {
            aceton_value_max += aceton.at<float>(i,j);
        }
    }
    aceton_value_max /= aceton.cols * aceton.rows;
    
    for (i = 0; i < aceton_min.rows; i++) {
        for (j = 0; j < aceton_min.cols; j++) {
            aceton_value_min += aceton_min.at<ushort>(i,j);
        }
    }
    aceton_value_min /= aceton_min.cols * aceton_min.rows;

    k = acenton_conc / (aceton_value_max - aceton_value_min); // УЖЕ с учетом 10000

    
    while (img_number < 501) {
        Mat PhotoLIF = imread(getNameAndFolder(img_number, input), CV_LOAD_IMAGE_ANYDEPTH);
        for (int i = 0; i < PhotoLIF.rows; i++) {
            for (int j = 0; j < PhotoLIF.cols; j++) {
                if ((PhotoLIF.at<ushort>(i,j) - background.at<ushort>(i,j)) < 0) {
                    PhotoLIF.at<ushort>(i,j) = 0;
                }
                else
                    PhotoLIF.at<ushort>(i,j) = PhotoLIF.at<ushort>(i,j) - background.at<ushort>(i,j);
            }
        }
        
        
        
        double temp = 0;
        
        for (int i = 0; i < PhotoLIF.rows; i++) {
            for (int j = 0; j < PhotoLIF.cols; j++) {
                
                   temp = static_cast<double>(PhotoLIF.at<ushort>(i,j)) / new_sheet.at<ushort>(i,j);
                   PhotoLIF.at<ushort>(i,j) = temp * 10000;
            }
        }
        
        for (int i = 0; i < PhotoLIF.rows; i++) {
            for (int j = 0; j < PhotoLIF.cols; j++) {
                PhotoLIF.at<ushort>(i,j) *= k * 10000000; // Концентрация будет * 10^-7
            }
        }
        
        
        output_image_name = getNameAndFolder(img_number, output);
        imwrite(output_image_name, PhotoLIF);
        img_number++;
    }

    return 0;
}


/*
 for (i = 1660; i < 1890; i++) {
 for (j = 0; j < new_sheet.cols; j++) {
 if (_new_sheet.at<float>(i,j) != 0)
 _new_sheet.at<float>(i,j) = static_cast<float>(PhotoLIF.at<ushort>(i,j)) /_new_sheet.at<float>(i,j);
 }
 }
 */
