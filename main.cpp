#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

typedef struct NeighborhoodStructure{
    int size;
    int* di;
    int* dj;
} NeighborhoodStructure;

typedef struct Kernel {
    float **values;
    int lengths;
    float meanValue;
} Kernel;

enum FilterType {
    minimum,
    median,
    maximum
};

float** initMatrix(int size) {
    float **matrix = (float**)calloc(size, sizeof(float*));

    for (int i = 0; i < size; i++) {
        matrix[i] = (float*)calloc(size, sizeof(float));
    }

    return matrix;
}

void printMatrix(int size, float **matrix) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << matrix[i][j] << " ";
        }

        cout << "\n";
    }
}

Mat applyOrderedFilter(Mat source, FilterType filter, int w) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));
    int middleW = w / 2;

    int loc = 0;

    for (int i = middleW; i < rows - middleW; i++) {
        for (int j = middleW; j < cols - middleW; j++) {
            vector<int> window;

            for (int k = i - middleW; k <= i + middleW; k++) {
                for (int l = j - middleW; l <= j + middleW; l++) {
                    window.push_back(source.at<uchar>(k, l));
                }
            }

            sort(window.begin(), window.end());

            if (filter == median) {
                loc = window.size() / 2;
            } else if (filter == maximum) {
                loc = window.size() - 1;
            }

            destination.at<uchar>(i, j) = window.at(loc);
        }
    }

    return destination;
}

Kernel computeKernel(int w) {
    int x0 = w / 2, y0 = w / 2;
    float sigma = 0.8;
    float pi = 3.14;
    float sum = 0.0;

    Kernel kernel;
    kernel.values = initMatrix(w);

    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            float x = (float)(pow(i - x0, 2) + (float)pow(j - y0, 2));
            float y = 2 * (float)pow(sigma, 2);
            float expVal = x / y;
            float a = exp(-expVal);
            float b = (2 * pi * (float)pow(sigma, 2));
            float val = a / b;

            sum += val;
            kernel.values[i][j] = val;
        }
    }

    kernel.lengths = w;
    kernel.meanValue = sum;

    return kernel;
}

Mat applyGaussianFilter(Mat source, Kernel kernel) {
    double t = (double)getTickCount();

    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));
    int border = kernel.lengths / 2;

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            float sum = 0;

            for (int k = 0; k < kernel.lengths; k++) {
                int newI = i + k - border;

                for (int l = 0; l < kernel.lengths; l++) {
                    int newJ = j + l - border;

                    float kernelValue = kernel.values[k][l];
                    float sourceValue = (float)source.at<uchar>(newI, newJ);

                    sum += kernelValue * sourceValue;
                }
            }

            destination.at<uchar>(i, j) = (uchar)max(min(sum / kernel.meanValue, 255.0f), 0.0f);
        }
    }

    t = ((double)getTickCount() - t) / getTickFrequency();
    printf("Time = %.3f [ms] to apply the Gaussian convolution with a kernel of size %d\n", t * 1000, kernel.lengths);

    return destination;
}

int main() {
    Mat balloonsSalty = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 10/PI-L10/balloons_Salt&Pepper.bmp",
                           IMREAD_GRAYSCALE);
    Mat balloonsGaussian = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 10/PI-L10/balloons_Gauss.bmp",
                          IMREAD_GRAYSCALE);
    imshow("Original Image Salt and Pepper", balloonsSalty);

    Mat orderFiltered = applyOrderedFilter(balloonsSalty, median, 7);
    imshow("Order Filtered Image", orderFiltered);

    Kernel kernel = computeKernel(7);
//    printMatrix(7, kernel.values);

    imshow("Original Image Gaussian Noise", balloonsGaussian);
    Mat gaussianFiltered = applyGaussianFilter(balloonsGaussian, kernel);
    imshow("Gaussian Filtered Image", gaussianFiltered);

    waitKey();

    return 0;
}
