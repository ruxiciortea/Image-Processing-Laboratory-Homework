//
// Created by Ruxandra Ciortea on 17.05.2023.
//

#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

typedef struct Kernel {
    float **values;
    int lengths;
    float meanValue;
} Kernel;

float** initMatrix(int rows, int cols);
void printMatrix(int rows, int cols, float **matrix);
int findRegion(double angle);
void findPointsBasedOnRegion(int region, Point &p1, Point &p2);
Kernel initKernel(vector<int> values, int size);
Mat applyThresholding(Mat source, int thresholdValue);
