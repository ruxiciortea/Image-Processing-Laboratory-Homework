#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

///Users/ruxiciortea/Desktop/IP/Lab 1/PI-L1/cell.bmp

Mat negative_image(Mat image) {
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    Mat outImage(image.rows, image.cols, CV_8UC1);

    for (int i = 0; i < image.rows; i++) {
        for (int j = 0; j < image.cols; j++) {
            outImage.at<uchar>(i,j) = 255 - image.at<uchar>(i,j);
        }
    }

    return outImage;

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****
}

Mat add_scalar(Mat image, int factor) {
    /*
     * Add a scalar to the entire image
     * Hint:
     *  Values may overshoot
     *  Checkout clone() method provided by OpenCV
     * Args:
     *  image: input grayscale image
     *  factor: value to add to each pixel of the image
     * Variables:
     *  rows: number of rows of the image
     *  cols: number of columns of the image
     * Return:
     *  result: transformed matrix
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    int rows = image.rows;
    int cols = image.cols;
    Mat result = image.clone();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uchar newVal = result.at<uchar>(i,j) + factor;

            if (newVal > 255) {
                result.at<uchar>(i, j) = 255;
            } else {
                result.at<uchar>(i, j) = newVal;
            }

            if (newVal < 0) {
                result.at<uchar>(i, j) = 0;
            } else {
                result.at<uchar>(i, j) = newVal;
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return result;
}

Mat mul_scalar(Mat image, float factor) {
    /*
     * Multiply a scalar with the entire image
     * Hint:
     *  Values may overshoot
     *  Checkout clone() method provided by OpenCV
     * Args:
     *  image: input grayscale image
     *  factor: value to multiply to each pixel of the image
     * Variables:
     *  rows: number of rows of the image
     *  cols: number of columns of the image
     * Return:
     *  result: transformed matrix
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    int rows = image.rows;
    int cols = image.cols;
    Mat result = image.clone();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uchar newVal = result.at<uchar>(i,j) * factor;

            if (newVal > 255) {
                result.at<uchar>(i, j) = 255;
            } else {
                result.at<uchar>(i, j) = newVal;
            }

            if (newVal < 0) {
                result.at<uchar>(i, j) = 0;
            } else {
                result.at<uchar>(i, j) = newVal;
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return result;

}

Mat draw_squares(int rows, int cols) {
    /*
     * Initialize a Mat object in order to create a square divided in four sub-squares which you are
     * going to color from top to bottom, left to right as: white, red, green, yellow
     * Hint:
     *  The channels are BGR not RGB
     *  value 0 means black 255 value means white (0 intensity to full intensity)
     *  You can initialize the image with White at the beginning
     * Args:
     *  rows: number of rows of the image
     *  cols: number of cols of the image
     * Variables:
     *  red: vector representing red color
     *  green: vector representing green color
     *  yellow: vector representing yellow color
     * Return:
     *  result: final matrix
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    Mat result(rows, cols, CV_8UC3);

    // color upper left corner
    Vec3b white(255, 255, 255);
    for (int i = 0; i < rows / 2 - 3; ++i) {
        for (int j = 0; j < cols / 2 - 3; ++j) {
            result.at<Vec3b>(i, j) = white;
        }
    }

    // color upper right corner
    Vec3b red(0, 0, 255);
    for (int i = 0; i < rows / 2 - 3; ++i) {
        for (int j = cols / 2 + 3; j < cols; ++j) {
            result.at<Vec3b>(i, j) = red;
        }
    }

    // color lower left corner
    Vec3b green(0, 255, 0);
    for (int i = rows / 2 + 3; i < rows; ++i) {
        for (int j = 0; j < cols / 2 - 3; ++j) {
            result.at<Vec3b>(i, j) = green;
        }
    }

    // color lower right corner
    Vec3b yellow(0, 255, 255);
    for (int i = rows / 2 + 3; i < rows; ++i) {
        for (int j = cols / 2 + 3; j < cols; ++j) {
            result.at<Vec3b>(i, j) = yellow;
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return result;
}

int main() {
    Mat image = imread("/Users/ruxiciortea/Desktop/IP/Lab 1/PI-L1/cell.bmp",
                       IMREAD_GRAYSCALE);

    Mat mul_image = mul_scalar(image, 0.7);

    imshow("Input image", image);
    imshow("Negative image", negative_image(image));
    imshow("Add image", add_scalar(image, 50));
    imshow("Mul image", mul_image);
    imshow("Drawing", draw_squares(256, 256));

    waitKey(0);
    return 0;
}
