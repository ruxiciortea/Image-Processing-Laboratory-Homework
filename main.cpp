#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

typedef struct Kernel {
    float **values;
    int lengths;
    float meanValue;
} Kernel;

typedef struct RowKernel {
    float *values;
    int length;
    float meanValue;
} RowKernel;

typedef struct ColumnKernel {
    float **values;
    int height;
    float meanValue;
} ColumnKernel;

enum FilterType {
    minimum,
    median,
    maximum
};

float** initMatrix(int rows, int cols) {
    float **matrix = (float**)calloc(rows, sizeof(float*));

    for (int i = 0; i < rows; i++) {
        matrix[i] = (float*)calloc(cols, sizeof(float));
    }

    return matrix;
}

void printMatrix(int rows, int cols, float **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << matrix[i][j] << " ";
        }

        cout << "\n";
    }

    cout << "\n";
}

void printArray(int elems, float *array) {
    for (int i = 0; i < elems; i++) {
        cout << array[i] << " ";
    }

    cout << "\n\n";
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
    float sigma = (float)w / 6;
    float pi = 3.14;
    float sum = 0.0;

    Kernel kernel;
    kernel.values = initMatrix(w, w);

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

RowKernel computeRowKernel(int w) {
    int x0 = w / 2;
    float sigma = (float)w / 6;
    float pi = 3.14;
    float sum = 0.0;

    RowKernel kernel;
    kernel.values = (float*)calloc(w, sizeof(float));

    for (int i = 0; i < w; i++) {
        float x = (float)pow(i - x0, 2);
        float y = 2 * (float)pow(sigma, 2);
        float expVal = x / y;
        float a = exp(-expVal);
        float b = sqrt(2 * pi * sigma);
        float val = a / b;

        sum += val;
        kernel.values[i] = val;
    }

    kernel.length = w;
    kernel.meanValue = sum;

    return kernel;
}

ColumnKernel computeColumnKernel(int w) {
    int y0 = w / 2;
    float sigma = (float)w / 6;
    float pi = 3.14;
    float sum = 0.0;

    ColumnKernel kernel;
    kernel.values = initMatrix(w, 1);

    for (int i = 0; i < w; i++) {
        float x = (float)pow(i - y0, 2);
        float y = 2 * (float)pow(sigma, 2);
        float expVal = x / y;
        float a = exp(-expVal);
        float b = sqrt(2 * pi * sigma);
        float val = a / b;

        sum += val;
        kernel.values[i][0] = val;
    }

    kernel.height = w;
    kernel.meanValue = sum;

    return kernel;
}

Mat applyGaussianFilter2D(Mat source, Kernel kernel) {
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
    printf("Time = %.3f [ms] to apply the Gaussian 2D convolution with a kernel of size %d\n", t * 1000, kernel.lengths);

    return destination;
}

Mat applyGaussianFilter1D(Mat source, RowKernel rowKernel, ColumnKernel columnKernel) {
    double t = (double)getTickCount();

    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));
    Mat aux = Mat(rows, cols, CV_8UC1, Scalar(0));

    int border = rowKernel.length / 2;

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            float sum = 0;

            for (int k = 0; k < columnKernel.height; k++) {
                int newJ = j + k - border;

                float kernelValue = columnKernel.values[k][0];
                float sourceValue = (float)source.at<uchar>(i, newJ);

                sum += kernelValue * sourceValue;
            }

            aux.at<uchar>(i, j) = (uchar)max(min(sum / columnKernel.meanValue, 255.0f), 0.0f);
        }
    }

    for (int i = border; i < cols - border; i++) {
        for (int j = border; j < rows - border; j++) {
            float sum = 0;

            for (int k = 0; k < rowKernel.length; k++) {
                int newI = i + k - border;

                float kernelValue = rowKernel.values[k];
                float sourceValue = (float)aux.at<uchar>(j, newI);

                sum += kernelValue * sourceValue;
            }

            destination.at<uchar>(j, i) = (uchar)max(min(sum / rowKernel.meanValue, 255.0f), 0.0f);
        }
    }

    t = ((double)getTickCount() - t) / getTickFrequency();
    printf("Time = %.3f [ms] to apply the Gaussian 1D convolution with a kernel of size %d\n", t * 1000, rowKernel.length);

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

    int kernel2DSize = 7;
    Kernel kernel = computeKernel(kernel2DSize);
    printMatrix(kernel2DSize, kernel2DSize, kernel.values);

    imshow("Original Image Gaussian Noise", balloonsGaussian);
    Mat gaussianFiltered2D = applyGaussianFilter2D(balloonsGaussian, kernel);
    imshow("Gaussian Filtered Image 2D", gaussianFiltered2D);

    int kernel1DSize = 7;
    RowKernel rowKernel = computeRowKernel(kernel1DSize);
    printArray(kernel1DSize, rowKernel.values);

    ColumnKernel columnKernel = computeColumnKernel(kernel1DSize);
    printMatrix(kernel1DSize, 1, columnKernel.values);

    Mat gaussianFiltered1D = applyGaussianFilter1D(balloonsGaussian, rowKernel, columnKernel);
    imshow("Gaussian Filtered Image 1D", gaussianFiltered1D);

    waitKey();

    return 0;
}
