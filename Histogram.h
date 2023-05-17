//
// Created by Ruxandra Ciortea on 17.05.2023.
//

#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

#define HISTOGRAM_SIZE 256

int* computeHistogram(Mat source);
void showHistogram(const string& name, int* hist, const int  hist_cols, const int hist_height);
