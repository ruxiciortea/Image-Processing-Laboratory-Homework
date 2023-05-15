#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

typedef struct NeighborhoodStructure{
    int size;
    int* di;
    int* dj;
} NeighborhoodStructure;

enum FilterType {
    minimum,
    median,
    maximum
};

Mat applyOrderedFilter(Mat source, FilterType filter, NeighborhoodStructure neighborhood) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));

    int loc = 0;

    for (int i = 1; i <= rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            vector<int> window;

            for (int k = 0; k < neighborhood.size; k++) {
                for (int l = 0; l < neighborhood.size; l++) {
                    window.push_back(source.at<uchar>(i + neighborhood.di[k], j + neighborhood.dj[l]));
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

int main() {
    Mat balloons = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 10/PI-L10/balloons_Salt&Pepper.bmp",
                           IMREAD_GRAYSCALE);
    imshow("Original Image", balloons);

    NeighborhoodStructure n8 = {
            8,
            (int[]){0,-1,-1,-1,0,1,1,1},
            (int[]){1,1,0,-1,-1,-1,0,1}
    };

    Mat orderFiltered = applyOrderedFilter(balloons, median, n8);
    imshow("Order Filtered Image", orderFiltered);

    waitKey();

    return 0;
}
