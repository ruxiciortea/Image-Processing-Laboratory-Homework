#include <iostream>
#include <opencv2/opencv.hpp>
#include "Histogram.h"
#include "Utils.h"

using namespace std;
using namespace cv;

#define STRONG_VALUE 255
#define WEAK_VALUE 128

Mat applyConvolution(Mat source, Kernel kernel, int border) {
    int rows = source.rows, cols = source.cols;
    Mat destination = source.clone();

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            float sum = 0;

            for (int k = 0; k < kernel.lengths; k++) {
                int newI = i + k - border - 1;

                for (int l = 0; l < kernel.lengths; l++) {
                    int newJ = j + l - border - 1;

                    float kernelValue = kernel.values[k][l];
                    float sourceValue = (float)source.at<uchar>(newI, newJ);

                    sum += kernelValue * sourceValue;
                }
            }

            destination.at<uchar>(i, j) = (uchar)max(min((float)sum / kernel.meanValue, 255.0f), 0.0f);
        }
    }

    return destination;
}

Mat computeMagnitudeFloat(Mat source, Kernel kernelX, Kernel kernelY) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_32FC1, Scalar(0));
    int border = 2;

    Mat filteredX = applyConvolution(source, kernelX, border);
    Mat filteredY = applyConvolution(source, kernelY, border);

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
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
            destination.at<uchar>(i, j) = (uchar)source.at<float>(i, j) / (4.0 * sqrt(2));
        }
    }

    return destination;
}

Mat computeDirection(Mat source, Kernel kernelX, Kernel kernelY) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_32FC1, Scalar(0));
    float pi = 3.14;
    int border = 2;

    Mat filteredX = applyConvolution(source, kernelX, border);
    Mat filteredY = applyConvolution(source, kernelY, border);

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            uchar xGrad = filteredX.at<uchar>(i, j);
            uchar yGrad = filteredY.at<uchar>(i, j);

            float angle = (float)(atan2(yGrad, xGrad) * (180.0 / pi));

            destination.at<float>(i, j) = angle;
        }
    }

    return destination;
}

Mat computeNonMaximaSuppression(Mat magnitudeMatrix, Mat directionMatrix) {
    int rows = magnitudeMatrix.rows, cols = magnitudeMatrix.cols;
    Mat destination = Mat(rows, cols, CV_32FC1, Scalar(0));
    int border = 2;

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            double direction = directionMatrix.at<uchar>(i, j);
            int region = findRegion(direction);

            Point p1, p2;
            findPointsBasedOnRegion(region, p1, p2);

            float m = magnitudeMatrix.at<float>(i, j);
            float m1 = magnitudeMatrix.at<float>(i + p1.x, j + p1.y);
            float m2 = magnitudeMatrix.at<float>(i + p2.x, j + p2.y);

            if (m > m1 && m > m2 && m > 0 && m1 > 0 && m2 > 0) {
                destination.at<float>(i, j) = magnitudeMatrix.at<float>(i, j);
            }
        }
    }

    return destination;
}

int computeAdaptiveThresholdValue(Mat source, int *magnitudeHistogramScaled) {
    int rows = source.rows, cols = source.cols;
    float p = 0.1;
    int border = 2;
    int numberNoEdges = (int)((1 - p) * (float)((rows - border) * (cols - border) - magnitudeHistogramScaled[0]));

    int adaptiveThreshold = 1;
    bool foundThreshold = false;
    int histSum = 0;

    while (adaptiveThreshold < HISTOGRAM_SIZE && !foundThreshold) {
        histSum += magnitudeHistogramScaled[adaptiveThreshold];

        if (histSum >= numberNoEdges) {
            foundThreshold = true;
        }

        adaptiveThreshold++;
    }

    return adaptiveThreshold;
}

Mat edgeLabeling(Mat source, int adaptiveThreshold) {
    float k = 0.4;
    int thresholdHigh = adaptiveThreshold;
    int thresholdLow = (int)(k * (float)thresholdHigh);

    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));
    int border = 1;

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            uchar pixelValue = (uchar)source.at<float>(i, j);

            if (pixelValue > thresholdHigh) {
                destination.at<uchar>(i, j) = STRONG_VALUE;
            } else if (pixelValue < thresholdHigh && pixelValue > thresholdLow) {
                destination.at<uchar>(i, j) = WEAK_VALUE;
            }
        }
    }

    return destination;
}

Mat extendStrongEdges(Mat source) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));
    Mat visitedPoints = Mat(rows, cols, CV_8UC1, Scalar(0));
    int border = 2;

    int neighborI[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int neighborJ[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int neighborSize = 8;

    queue<Point> queue;

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            if (source.at<uchar>(i, j) == STRONG_VALUE && visitedPoints.at<uchar>(i, j) != STRONG_VALUE) {
                queue.push(Point(i, j));

                visitedPoints.at<uchar>(i, j) = STRONG_VALUE;
                destination.at<uchar>(i, j) = STRONG_VALUE;
            }

            while (!queue.empty()) {
                Point p = queue.front();
                queue.pop();

                for (int k = 0; k < neighborSize; k++) {
                    int newI = p.x+ neighborI[k];

                    for (int l = 0; l < neighborSize; l++) {
                        int newJ = p.y + neighborJ[l];

                        if (source.at<uchar>(newI, newJ) == WEAK_VALUE && visitedPoints.at<uchar>(newI, newJ) != 255) {
                            queue.push(Point(newI, newJ));

                            visitedPoints.at<uchar>(newI, newJ) = STRONG_VALUE;
                            destination.at<uchar>(newI, newJ) = STRONG_VALUE;
                        }
                    }
                }
            }
        }
    }

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            if (source.at<uchar>(i, j) == WEAK_VALUE) {
                destination.at<uchar>(i, j) = 0;
            }
        }
    }

    return destination;
}

int main() {
    Mat saturn = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 11/PI-L11/saturn.bmp",
                               IMREAD_GRAYSCALE);

    int kernelSize3 = 3;

    Kernel gaussianKernel = initKernel({1, 2, 1, 2, 4, 2, 1, 2, 1}, kernelSize3);
    Kernel sobelX = initKernel({-1, 0, 1, -2, 0, 2, -1, 0, 1}, kernelSize3);
    Kernel sobelY = initKernel({1, 2, 1, 0, 0, 0, -1, -2, -1}, kernelSize3);

    imshow("Original Image", saturn);
    Mat gaussianFiltered2D = applyConvolution(saturn, gaussianKernel, 1);
    imshow("Gaussian Filtered Image", gaussianFiltered2D);

    Mat sobelXFiltered = applyConvolution(gaussianFiltered2D, sobelX, 2);
    Mat sobelYFiltered = applyConvolution(gaussianFiltered2D, sobelY, 2);

    Mat floatMagnitudeMatrix = computeMagnitudeFloat(gaussianFiltered2D, sobelX, sobelY);
    Mat normalizedMagnitudeMatrix = normalizeMagnitude(floatMagnitudeMatrix);
    imshow("Magnitude values", normalizedMagnitudeMatrix);

    Mat directionMatrix = computeDirection(gaussianFiltered2D, sobelX, sobelY);

    Mat nonMaximaSuppression = computeNonMaximaSuppression(floatMagnitudeMatrix, directionMatrix);
    Mat normalizedNonMaximaSuppression;
    normalize(nonMaximaSuppression, normalizedNonMaximaSuppression, 0, 1, NORM_MINMAX);
    imshow("Normalized Non Maxima Superposition", normalizedNonMaximaSuppression);

    int *magnitudeHistogramScales = computeHistogram(normalizedMagnitudeMatrix);
//    showHistogram("Histogram", magnitudeHistogramScales, HISTOGRAM_SIZE, 100);
    int adaptiveThreshold = computeAdaptiveThresholdValue(nonMaximaSuppression, magnitudeHistogramScales);

    Mat thresholding = edgeLabeling(nonMaximaSuppression, adaptiveThreshold);
    imshow("After Thresholding", thresholding);

    Mat edges = extendStrongEdges(thresholding);
    imshow("Just edges", edges);

    waitKey();

    return 0;
}
