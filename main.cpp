#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

#define HISTOGRAM_SIZE 256

int* compute_histogram(Mat source, int histogram_bins) {

    int rows = source.rows;
    int cols = source.cols;
    int* histogram = new int[HISTOGRAM_SIZE];

    for (int i = 0; i < HISTOGRAM_SIZE; ++i) {
        histogram[i] = 0;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            unsigned char pixel = source.at<uchar>(i, j);
            histogram[pixel]++;
        }
    }

    return histogram;

}

void showHistogram(const string& name, int* hist, const int  hist_cols, const int hist_height) {

    // constructs a white image
    Mat imgHist(hist_height, hist_cols, CV_8UC3, CV_RGB(255, 255, 255));

    //computes histogram maximum
    int max_hist = 0;

    for (int i = 0; i < hist_cols; i++) {
        if (hist[i] > max_hist) {
            max_hist = hist[i];
        }
    }

    double scale = (double)hist_height / max_hist;
    int baseline = hist_height - 1;

    for (int x = 0; x < hist_cols; x++) {
        Point p1 = Point(x, baseline);
        Point p2 = Point(x, baseline - cvRound(hist[x] * scale)); line(imgHist, p1, p2, CV_RGB(255, 0, 255)); // histogram bins
    }

    imshow(name, imgHist);

}

float compute_mean_intensity(Mat source, int* histogram) {

    int rows = source.rows, cols = source.cols;
    int M = rows * cols;
    int sum = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            sum += source.at<uchar>(i, j);
        }
    }

    return (float)sum / M;

}

float compute_standard_deviation(Mat source, int* histogram) {

    int rows = source.rows, cols = source.cols;
    int M = rows * cols;
    int sum = 0;
    float mean = compute_mean_intensity(source, histogram);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float elem = source.at<uchar>(i, j) - mean;
            sum += elem * elem;
        }
    }

    return sqrt(float(sum) / M);

}

int compute_cumulative_histogram(int* histogram, int start_level, int end_level) {

    int rez = 0;

    for (int i = start_level; i < end_level; i++) {
        rez += histogram[i];
    }

    return rez;

}

int histogram_min(int* histogram) {

    for (int i = 0; i < HISTOGRAM_SIZE; i++) {
        if (histogram[i] != 0) {
            return i;
        }
    }

    return 0;

}

int histogram_max(int* histogram) {

    for (int i = HISTOGRAM_SIZE - 1; i >= 0; i--) {
        if (histogram[i] != 0) {
            return i;
        }
    }

    return HISTOGRAM_SIZE;

}

int compute_global_threshold(Mat source, int* histogram) {

    int min = histogram_min(histogram);
    int max = histogram_max(histogram);
    int T = (min + max) / 2;

    int N1 = compute_cumulative_histogram(histogram, min, T);
    int N2 = compute_cumulative_histogram(histogram, T, max);

    int sumG1 = 0;
    int sumG2 = 0;

    for (int i = min; i <= T; i++) {
        sumG1 += i * histogram[i];
    }

    for (int i = T + 1; i <= max; i++) {
        sumG2 += i * histogram[i];
    }

    float meanG1 = (float)sumG1 / N1;
    float meanG2 = (float)sumG2 / N2;

    return (meanG1 + meanG2) / 2;

}

Mat grayscale_2_binary(Mat source, int threshold){

    int rows = source.rows, cols = source.cols;
    Mat binary(rows, cols, CV_8UC1);;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            unsigned char pixel = source.at<unsigned char>(i,j);

            if (pixel <= threshold) {
                binary.at<uchar>(i, j) = 0;
            } else {
                binary.at<uchar>(i, j) = 255;
            }
        }
    }

    return binary;

}

Mat histogram_slide(Mat source, int offset, int min, int max) {

    int rows = source.rows, cols = source.cols;
    Mat dst(rows, cols, CV_8UC1, Scalar(0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int new_value = source.at<uchar>(i, j) + offset;

            if (new_value < min) {
                dst.at<uchar>(i, j) = min;
            } else if (new_value > max) {
                dst.at<uchar>(i, j) = max;
            } else {
                dst.at<uchar>(i, j) = new_value;
            }
        }
    }

    return dst;

}

Mat histogram_shrinking_stretching(Mat source, int g_out_min, int g_out_max) {

    int rows = source.rows, cols = source.cols;
    Mat dst(rows, cols, CV_8UC1, Scalar(0));

    int g_in_min = 256;
    int g_in_max = -1;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uchar current_pixel = source.at<uchar>(i, j);

            if (current_pixel < g_in_min) {
                g_in_min = current_pixel;
            }

            if (current_pixel > g_in_max) {
                g_in_max = current_pixel;
            }
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float fraction = (float)(g_out_max - g_out_min) / (g_in_max - g_in_min);
            int new_value = g_out_min + (source.at<uchar>(i, j) - g_in_min) * fraction;

            if (new_value < g_out_min) {
                dst.at<uchar>(i, j) = g_out_min;
            } else if (new_value > g_out_max) {
                dst.at<uchar>(i, j) = g_out_max;
            } else {
                dst.at<uchar>(i, j) = new_value;
            }
        }
    }

    return dst;

}

int main() {

    Mat source = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 8/PI-L8/balloons.bmp",IMREAD_GRAYSCALE);
    Mat nature = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 8/PI-L8/Hawkes_Bay_NZ.bmp",IMREAD_GRAYSCALE);
    Mat wheel = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 8/PI-L8/wheel.bmp",IMREAD_GRAYSCALE);

    imshow("Original Image", source);

    int* histogram_balloon = compute_histogram(source, HISTOGRAM_SIZE);
    showHistogram("Histogram", histogram_balloon, HISTOGRAM_SIZE, 100);

    float mean = compute_mean_intensity(source, histogram_balloon);
    float deviation = compute_standard_deviation(source, histogram_balloon);
    int cumulativeHistogram = compute_cumulative_histogram(histogram_balloon, 0, 120);

    printf("Mean intensity: %f \n", mean);
    printf("Standard deviation: %f \n", deviation);
    printf("Cumulative histogram: %d \n", cumulativeHistogram);

//    int threshold = compute_global_threshold(source, histogram_balloon);
//    Mat global_threshold_binary_image = grayscale_2_binary(source, threshold);
//    imshow("After Thresholding", global_threshold_binary_image);

//    Mat balloons_brighter = histogram_slide(source, 50, 0, 255);
//    Mat balloons_darker = histogram_slide(source, -50, 0, 255);
//
//    imshow("Balloons Brighter", balloons_brighter);
//    imshow("Balloons Darker", balloons_darker);

    Mat nature_stretching = histogram_shrinking_stretching(nature, 10, 255);
    Mat wheel_shrinking = histogram_shrinking_stretching(wheel, 50, 150);

    imshow("Nature Original", nature);
    imshow("Nature Stretching", nature_stretching);
    imshow("Wheel Original", wheel);
    imshow("Wheel Shrinking", wheel_shrinking);

    waitKey();

    return 0;

}
