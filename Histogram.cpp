#include "Histogram.h"

int* computeHistogram(Mat source) {
    int rows = source.rows;
    int cols = source.cols;
    int* histogram = new int[HISTOGRAM_SIZE];

    for (int i = 0; i < HISTOGRAM_SIZE; ++i) {
        histogram[i] = 0;
    }

    for (int i = 2; i < rows - 2; i++) {
        for (int j = 2; j < cols - 2; j++) {
            unsigned char pixel = source.at<uchar>(i, j);
            histogram[pixel]++;
        }
    }

    return histogram;
}

void showHistogram(const string& name, int* hist, const int  hist_cols, const int hist_height) {
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