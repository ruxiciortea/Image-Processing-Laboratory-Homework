#include "Histogram.h"

int* computeHistogram(Mat source) {
    int rows = source.rows;
    int cols = source.cols;
    int* histogram = (int*)calloc(HISTOGRAM_SIZE, sizeof(int));
    int border = 2;

    for (int i = border; i < rows - border; i++) {
        for (int j = border; j < cols - border; j++) {
            uchar pixel = source.at<uchar>(i, j);
            histogram[pixel]++;
        }
    }

    return histogram;
}

void showHistogram(const string& name, int* hist, const int  histCols, const int histHeight) {
    Mat imgHist(histHeight, histCols, CV_8UC3, CV_RGB(255, 255, 255));

    //computes histogram maximum
    int max_hist = 0;

    for (int i = 0; i < histCols; i++) {
        if (hist[i] > max_hist) {
            max_hist = hist[i];
        }
    }

    double scale = (double)histHeight / max_hist;
    int baseline = histHeight - 1;

    for (int x = 0; x < histCols; x++) {
        Point p1 = Point(x, baseline);
        Point p2 = Point(x, baseline - cvRound(hist[x] * scale)); line(imgHist, p1, p2, CV_RGB(255, 0, 255)); // histogram bins
    }

    imshow(name, imgHist);
}