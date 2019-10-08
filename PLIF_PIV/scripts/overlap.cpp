#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <numeric>

using namespace std;
using namespace cv;

const int total_images = 5459;

namespace PIV_coord {
    const double x_max(11.51885);
    const double x_min(-25.95448);
    const double y_max(20.53553);
    const double y_min(-17.28168);
    const double y_max_488(2.953842);
    const double pix_x(0.043371); // Сколько мм в 1 пикселе ПИВ. взял расстояние между двумя скоростями в PIV по х и поделил на 8
    const double pix_y(0.041466); // same только по у
    const double pix_in_LIF_x(1.15746); // PIVшный пиксель в LIV писелях по x
    const double pix_in_LIF_y(1.10662); // PIVшный пиксель в LIV писелях по у
}

namespace target {
    const double av_dot_dif(26.6875);
    const int x_center(705);
    const int y_center(555);
}

enum {lif_x_max, lif_x_min, lif_y_max, lif_y_min};

void createFolders() {
    boost::filesystem::create_directories("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_overlaped_noise_txt/");
    boost::filesystem::create_directories("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/Series_06_09_17_calibration/S1_Re5000_43A_T800_2D_25C_p1_C001H001S0001/temperature_overlaped_noise");
}

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

void getPixelValuePIVGrid(Mat& imagePhoto, vector<vector<vector<double>>> new_coordinate_matrix, vector<int> LIF_corners, string name) {
    int _j = LIF_corners[lif_y_max], _i = LIF_corners[lif_x_min];
    ofstream fout(name);
    
    fout << "x[mm]" << " " << "y[mm]" << " " << "T[Celsius]" << endl;
    
    for (double i = PIV_coord::x_min; i < PIV_coord::x_max; i += 8 * PIV_coord::pix_x) {
        for (double j = PIV_coord::y_max; j > PIV_coord::y_max_488; j -= 8 * PIV_coord::pix_y) {
            
            while(abs(new_coordinate_matrix[_i][0][0] - i) > PIV_coord::pix_x) {
                _i++;
            }
            
            while(abs(new_coordinate_matrix[0][_j][1] - j) > PIV_coord::pix_y) {
                _j++;
            }
            
            fout << i << " " << j << " " << imagePhoto.at<ushort>(_j - LIF_corners[lif_y_max], _i - LIF_corners[lif_x_min]) << endl;
            _i = LIF_corners[lif_x_min];
            _j = LIF_corners[lif_y_max];
        }
    }
}

void newCoordinateMatrix(vector<vector<vector<double>>> &coordinate_matrix, double pix_mm, vector<Vec3f> center) {
    pix_mm = 1 / pix_mm; //1 pix in mm
    ofstream output("/media/itpuser/IMP_JET_PROC_STOR1/Stereo/sayanLIF/matrix.txt");
    
    coordinate_matrix[center[0][0]][center[0][1]][0] = 0;
    coordinate_matrix[center[0][0]][center[0][1]][1] = 0;
    
    int j = 0;
    for (int x = center[0][0]; x > -1; x--) {
        int i = 0;
        for (int y = center[0][1]; y > -1; y--) {
            coordinate_matrix[x][y][1] = pix_mm * i;
            coordinate_matrix[x][y][0] = -pix_mm * j;
            i++;
        }
        j++;
    }
    
    int k = 0;
    for (int x = center[0][0]; x < 1024; x++) {
        int i = 0;
        for (int y = center[0][1]; y > -1; y--) {
            coordinate_matrix[x][y][1] = pix_mm * i;
            coordinate_matrix[x][y][0] = pix_mm * k;
            i++;
        }
        k++;
    }
    
    int l = 0;
    for (int x = center[0][0]; x < 1024; x++) {
        int i = 0;
        for (int y = center[0][1]; y < 1024; y++) {
            coordinate_matrix[x][y][1] = -pix_mm * i;
            coordinate_matrix[x][y][0] = pix_mm * l;
            i++;
        }
        l++;
    }
    
    int m = 0;
    for (int x = center[0][0]; x > -1; x--) {
        int i = 0;
        for (int y = center[0][1]; y < 1024; y++) {
            coordinate_matrix[x][y][1] = -pix_mm * i;
            coordinate_matrix[x][y][0] = -pix_mm * m;
            i++;
        }
        m++;
    }
    
    for (int _i = 0; _i < 1024; _i++) {
        for (int _j = 0; _j < 1024; _j++) {
            output << coordinate_matrix[_j][_i][0] << "," << coordinate_matrix[_j][_i][1] << " ";
        }
        output << endl;
    }
}

void findUtmostDots(vector<vector<vector<double>>> coordinate_matrix, vector<int> &LIF_corners, double pix_mm) {
    pix_mm = 1 / pix_mm; //1 pix in mm
    
    for (int _j = 0; _j < 1024; _j++) {
        if (abs(coordinate_matrix[_j][0][0] - PIV_coord::x_max) < pix_mm / 2) // сравниваем мм
            LIF_corners[lif_x_max] = _j + round(4 * PIV_coord::pix_in_LIF_x) - 1; //номер х в пиксельной матрице, 16 это из усреднения PIV по квадрату 32х32
        if (abs(coordinate_matrix[_j][0][0] - PIV_coord::x_min) < pix_mm / 2)
            LIF_corners[lif_x_min] = _j - round(4 * PIV_coord::pix_in_LIF_x) + 1;
    }
    
    for (int _j = 0; _j < 1024; _j++) {
        if (abs(coordinate_matrix[0][_j][1] - PIV_coord::y_max) < pix_mm / 2)
            LIF_corners[lif_y_max] = _j - round(4 * PIV_coord::pix_in_LIF_y); //номер y_мин в пиксельной матрице
        if (abs(coordinate_matrix[0][_j][1] - PIV_coord::y_max_488) < pix_mm / 2)
            LIF_corners[lif_y_min] = _j + round(4 * PIV_coord::pix_in_LIF_y); //номер y_мин в пиксельной матрице
    }
    
    //LIF_corners[lif_y_max] = 488; //тк лиф изображение у меня обрезано и имеет размер 1024х488
}

void averageOverlapLifImages(vector<int> LIF_corners, vector<vector<vector<double>>> new_coordinate_matrix) {
    int img_number = 1;
    bool input = 1, output = 0;
    string img_name;
    string output_image_name;
    
    Rect ROI(LIF_corners[lif_x_min], LIF_corners[lif_y_max], LIF_corners[lif_x_max] - LIF_corners[lif_x_min], LIF_corners[lif_y_min] - LIF_corners[lif_y_max]);
    
    while (img_number <= total_images) {
        Mat temperaturePhotoLIF = imread(getNameAndFolder(img_number, input), CV_LOAD_IMAGE_ANYDEPTH);
        temperaturePhotoLIF = temperaturePhotoLIF(ROI);
        double average_rect = 0;
        int m = 9, n = 9, _i = 1, _j = 1, _o = 0, __o = 0;
        
        for (int i = 0; i < temperaturePhotoLIF.rows; i += m) {
            if ((i + m) - PIV_coord::pix_in_LIF_y * 8 * _i > PIV_coord::pix_in_LIF_y) {
                m = 8;
            }
            else {
                m = 9;
            }
            for (int j = 0; j < temperaturePhotoLIF.cols; j += n) {
                
                if (PIV_coord::pix_in_LIF_x * 8 * _j - (j + n) > PIV_coord::pix_in_LIF_x) {
                    n = 10;
                }
                else {
                    n = 9;
                }
                
                for (int k = 0; k < m; k++) {
                    for (int l = 0; l < n; l++) {
                        average_rect += temperaturePhotoLIF.at<ushort>(i + k, j + l);
                    }
                }
                
                average_rect /= m * n;
                
                for (int k = 0; k < m; k++) {
                    for (int l = 0; l < n; l++)
                        temperaturePhotoLIF.at<ushort>(i + k, j + l) = average_rect;
                }
                average_rect = 0;
                _j++;
            }
            _j = 1;
            _i++;
        }
        
        // осреднеине сетки 16х16
        average_rect = 0;
        m = 9;
        n = 9;
        _i = 1;
        _j = 1;
        int temp_n = 0, temp_m = 0;
        
        for (int i = 0; i < temperaturePhotoLIF.rows; i += m) {
            
            average_rect = 0;
    
            if ((i + m) - PIV_coord::pix_in_LIF_y * 8 * _i > PIV_coord::pix_in_LIF_y) {
                m = 8;
            }
            else {
                m = 9;
            }
            temp_m = m;
            
            for (int j = 0; j < temperaturePhotoLIF.cols; j += n) {
                
                if (PIV_coord::pix_in_LIF_x * 8 * _j - (j + n) > PIV_coord::pix_in_LIF_x) {
                    n = 10;
                }
                else {
                    n = 9;
                }
                temp_n = n;
                
                for (int k = 0; k < 9 + m; k += temp_m) {
                    temp_n = n;
                    for (int l = 0; l < 9 + n; l += temp_n) {
                        average_rect += temperaturePhotoLIF.at<ushort>(i + k, j + l);
                        if (PIV_coord::pix_in_LIF_x * 8 * (_j + 1) - (j + 2 * n) > PIV_coord::pix_in_LIF_x) {
                            temp_n = 10;
                        }
                        else {
                            temp_n = n;
                        }
                        if (PIV_coord::pix_in_LIF_y * 8 * (_i + 1) - (i + 2 * m) > PIV_coord::pix_in_LIF_y) {
                            temp_m = 8;
                        }
                        else {
                            temp_m = 9;
                        }
                    }
                }
                average_rect /= 4;
                
                if (i >= 470) { //Проверка на квадраты снизу (осреднение с верхними)
                    average_rect = 0;
                    for (int k = 0; k < 9 + temperaturePhotoLIF.rows - i; k += temp_m) {
                        temp_n = n;
                        for (int l = 0; l < 9 + n; l += temp_n) {
                            average_rect += temperaturePhotoLIF.at<ushort>(i - 9 + k, j + l);
                        }
                    }
                    for (int k = 0; k < 8; k++) {
                        for (int l = 0; l < 9 + n; l++)
                        temperaturePhotoLIF.at<ushort>(i + k, j + l) = average_rect / 4;
                    }
                    _j++;
                    continue;
                }
                
                if (j >= 999)  { //Проверка на квадраты справа (осреднение с левыми)
                    average_rect = 0;
                    for (int k = 0; k < 9 + m; k += temp_m) {
                        temp_n = n;
                        for (int l = 0; l < 9 + temperaturePhotoLIF.cols - j; l += temp_n) {
                            average_rect += temperaturePhotoLIF.at<ushort>(i + k, j - 9 + l);
                        }
                    }
                    
                    for (int k = 0; k < 9 + m; k++) {
                        for (int l = 0; l < 10; l++)
                        temperaturePhotoLIF.at<ushort>(i + k, j - n + l) = average_rect / 4;
                    }
                    _j++;
                    continue;
                }
                
                for (int k = 0; k < 9 + m; k++) {
                    for (int l = 0; l < 9 + n; l++)
                        temperaturePhotoLIF.at<ushort>(i + k, j + l) = average_rect;
                }
                
                average_rect = 0;
                _j++;
            }
            _j = 1;
            _i++;
        }
        // 14.02 добавил ROI
        // temperaturePhotoLIF = temperaturePhotoLIF(ROI);
        
        imwrite(getNameAndFolder(img_number, output), temperaturePhotoLIF);
        getPixelValuePIVGrid(temperaturePhotoLIF, new_coordinate_matrix, LIF_corners, getNameAndFolder(img_number, 3));
        if (img_number % 20)
        cout << (img_number * 1.0 / total_images) * 100.0 << "%.." << endl;
        img_number++;
    }
}

int main() {
    createFolders();
    double average_dot_difference = target::av_dot_dif;
    vector<Vec3f> circles2;
    circles2.push_back(0);
    circles2[0][0] = target::x_center;
    circles2[0][1] = target::y_center;
    
    vector<vector<vector<double>>> new_coordinate_matrix(1024, vector<vector<double>>(1024, vector<double>(2, 0)));
    vector<int> LIF_corners {1023, 0, 1023, 0};
    
    newCoordinateMatrix(new_coordinate_matrix, average_dot_difference, circles2);
    
    findUtmostDots(new_coordinate_matrix, LIF_corners, average_dot_difference); //координаты для обрезки LIF изображения
    averageOverlapLifImages(LIF_corners, new_coordinate_matrix);
    
    return 0;
}

//ssh -oKexAlgorithms=+diffie-hellman-group1-sha1 -X itpuser@192.168.1.130
