#include <iostream>
#include <opencv2/opencv.hpp>
#include "Histogram.h"
#include "Utils.h"

using namespace std;
using namespace cv;

Kernel computeGaussianKernel(int w) {
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

Mat applyConvolution(Mat source, Kernel kernel) {
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

    return destination;
}

Mat computeMagnitudeFloat(Mat source, Kernel kernelX, Kernel kernelY) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_32FC1, Scalar(0));

    Mat filteredX = applyConvolution(source, kernelX);
    Mat filteredY = applyConvolution(source, kernelY);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uchar xGrad = filteredX.at<uchar>(i, j);
            uchar yGrad = filteredY.at<uchar>(i, j);

            destination.at<float>(i, j) = abs(sqrt(pow(xGrad, 2) + pow(yGrad, 2)));
        }
    }

    return destination;
}

Mat normalizeMagnitude(Mat source) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            destination.at<uchar>(i, j) = source.at<float>(i, j) / (4 * sqrt(2));
        }
    }

    return destination;
}

Mat computeDirection(Mat source, Kernel kernelX, Kernel kernelY) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_32FC1, Scalar(0));
    float pi = 3.14;

    Mat filteredX = applyConvolution(source, kernelX);
    Mat filteredY = applyConvolution(source, kernelY);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uchar xGrad = filteredX.at<uchar>(i, j);
            uchar yGrad = filteredY.at<uchar>(i, j);

            double angle = atan((float)yGrad/xGrad) * 180 / pi;

            if (kernelX.lengths == 2 || kernelY.lengths == 2) {
                angle += 135;
            }

            destination.at<float>(i, j) = angle;
        }
    }

    return destination;
}

Mat nonMaximaSuperposition(Mat magnitudeMatrix, Mat directionMatrix) {
    int rows = magnitudeMatrix.rows, cols = magnitudeMatrix.cols;
    Mat destination = Mat(rows, cols, CV_32FC1, Scalar(0));

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            double direction = directionMatrix.at<uchar>(i, j);
            int region = findRegion(direction);

            Point p1, p2;
            findPointsBasedOnRegion(region, p1, p2);

            if (magnitudeMatrix.at<uchar>(i, j) > magnitudeMatrix.at<uchar>(i + p1.x, j + p1.y)
                    && magnitudeMatrix.at<uchar>(i, j) > magnitudeMatrix.at<uchar>(i + p2.x, j + p2.y)) {
                destination.at<float>(i, j) = magnitudeMatrix.at<float>(i, j);
            } else {
                destination.at<float>(i, j) = 0.0f;
            }
        }
    }

    return destination;
}

int computeAdaptiveThresholdValue(Mat source, int *magnitudeHistogramScales) {
    int rows = source.rows, cols = source.cols;
    float p = 0.1;
    int numberNoEdges = (1 - p) * (rows - 2) * (cols - 2) * magnitudeHistogramScales[0];

    int adaptiveThreshold = 1;
    bool foundThreshold = false;
    int histSum = 0;

    while (adaptiveThreshold < HISTOGRAM_SIZE && !foundThreshold) {
        histSum += magnitudeHistogramScales[adaptiveThreshold];

        if (histSum >= numberNoEdges) {
            foundThreshold = true;
        }

        adaptiveThreshold++;
    }

    return adaptiveThreshold;
}

Mat applyThresholding(Mat source, int thresholdValue){
    int rows = source.rows, cols = source.cols;
    Mat destination(rows, cols, CV_8UC1);;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            unsigned char pixel = source.at<unsigned char>(i,j);

            if (pixel <= thresholdValue) {
                destination.at<uchar>(i, j) = 0;
            } else {
                destination.at<uchar>(i, j) = pixel;
            }
        }
    }

    return destination;
}

int main() {
    Mat saturn = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 11/PI-L11/saturn.bmp",
                               IMREAD_GRAYSCALE);

    int gaussianKernelSize = 3;
    int kernelSize3 = 3;
    int kernelSize2 = 2;
    Kernel gaussianKernel = computeGaussianKernel(gaussianKernelSize);

//    Kernel prewittX = initKernel({-1, 0, 1, -1, 0, 1, -1, 0, 1}, kernelSize3);
//    Kernel prewittY = initKernel({1, 1, 1, 0, 0, 0, -1, -1, -1}, kernelSize3);
    Kernel sobelX = initKernel({-1, 0, 1, -2, 0, 2, -1, 0, 1}, kernelSize3);
    Kernel sobelY = initKernel({1, 2, 1, 0, 0, 0, -1, -2, -1}, kernelSize3);
//    Kernel crossX = initKernel({1, 0, 0, -1}, kernelSize2);
//    Kernel crossY = initKernel({0, -1, 1, 0}, kernelSize2);

    printMatrix(3, 3, sobelX.values);
    printMatrix(3, 3, sobelY.values);

    imshow("Original Image", saturn);
    Mat gaussianFiltered2D = applyConvolution(saturn, gaussianKernel);
//    imshow("Gaussian Filtered Image", gaussianFiltered2D);

//    Mat prewittXFiltered = applyConvolution(gaussianFiltered2D, prewittX);
//    Mat prewittYFiltered = applyConvolution(prewittXFiltered, prewittY);
//    imshow("Prewitt Filtered Image", prewittYFiltered);

    Mat sobelXFiltered = applyConvolution(gaussianFiltered2D, sobelX);
    Mat sobelYFiltered = applyConvolution(gaussianFiltered2D, sobelY);
//    imshow("X", sobelXFiltered);
//    imshow("Y", sobelYFiltered);

//    Mat crossXFiltered = applyConvolution(gaussianFiltered2D, crossX);
//    Mat crossYFiltered = applyConvolution(crossXFiltered, crossY);
//    imshow("Cross Filtered Image", crossYFiltered);

    Mat floatMagnitudeMatrix = computeMagnitudeFloat(gaussianFiltered2D, sobelX, sobelY);
    Mat normalizedMagnitudeMatrix = normalizeMagnitude(floatMagnitudeMatrix);
    imshow("Magnitude values", normalizedMagnitudeMatrix);

    Mat directionMatrix = computeDirection(gaussianFiltered2D, sobelX, sobelY);
//    imshow("Direction values", directionMatrix);

    Mat nonMaximaSuperpositioned = nonMaximaSuperposition(floatMagnitudeMatrix, directionMatrix);
    Mat normalizedNonMaximaSuperpositioned;
    normalize(nonMaximaSuperpositioned, normalizedNonMaximaSuperpositioned, 0, 1, NORM_MINMAX);
    imshow("test", normalizedNonMaximaSuperpositioned);

    int *magnitudeHistogramScales = computeHistogram(normalizedMagnitudeMatrix);
//    showHistogram("Histogram", magnitudeHistogramScales, HISTOGRAM_SIZE, 100);
    int adaptiveThreshold = computeAdaptiveThresholdValue(nonMaximaSuperpositioned, magnitudeHistogramScales);

    Mat thresholding = applyThresholding(nonMaximaSuperpositioned, adaptiveThreshold);
    imshow("x", thresholding);

    waitKey();

    return 0;
}
