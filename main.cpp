#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

typedef struct kernel {
    int values[3][3];
    int lengths;
    int mean_value;
};

Mat filter(Mat source, kernel kernel) {
    int rows = source.rows, cols = source.cols;
    Mat destination = Mat(rows, cols, CV_8UC1, Scalar(0));

    for (int i = 1; i <= rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            int sum = 0;

            for (int k = 0; k < kernel.lengths; k++) {
                int new_i = i + k;

                for (int l = 0; l < kernel.lengths; l++) {
                    int new_j = j + k;

                    sum += kernel.values[k][l] * source.at<uchar>(new_i, new_j);
                }
            }

            destination.at<uchar>(i, j) = max(min(sum / kernel.mean_value, 255), 0);
        }
    }

    return destination;
}

void centering_transform(Mat source) {
    int rows = source.rows, cols = source.cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            source.at<float>(i, j) = ((i + j) & 1) ? - source.at<float>(i, j) : source.at<float>(i, j);
        }
    }
}

Mat generic_frequency_domain_filter(Mat source) {
    // convert the image to float image
    Mat source_float;
    source.convertTo(source_float, CV_32FC1);

    // apply the centering transformation
    centering_transform(source_float);

    // perform forward transform with complex image output
    Mat fourier;
    dft(source_float, fourier, DFT_COMPLEX_OUTPUT);

    // split real and imaginary channels
    Mat channels[] = {
            Mat::zeros(source.size(), CV_32F),
            Mat::zeros(source.size(), CV_32F)
    };
    split(fourier, channels); // channels[0] = Re(DFT(I)), channels[1] = Im(DFT(I))

    //calculate magnitude and phase in floating point images mag and phi
    Mat mag, phi;
    magnitude(channels[0], channels[1], mag);
    phase(channels[0], channels[1], phi);

    //display the phase and magnitude images here
    imshow("Magnitude", mag);
    imshow("Phase", phi);

    //insert filtering operations on Fourier coefficients here
    // ......

    //store in real part in channels[0] and imaginary part in channels[1]
    // ......

    //perform inverse transform and put results in dstf
    Mat dst, dstf;
    merge(channels, 2, fourier);
    dft(fourier, dstf, DFT_INVERSE | DFT_REAL_OUTPUT | DFT_SCALE);

    //inverse centering transformation
    centering_transform(dstf);

    //normalize the result and put in the destination image
    normalize(dstf, dst, 0, 255,  NORM_MINMAX, CV_8UC1);
    //Note: normalizing distorts the resut while enhancing the image display in the range [0,255].
    //For exact results (see Practical work 3) the normalization should be replaced with convertion:
    //dstf.convertTo(dst, CV_8UC1);

    return dst;
}

int main() {
    Mat cameraman = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 9/PI-L9/cameraman.bmp",
                        IMREAD_GRAYSCALE);

    imshow("Original Image", cameraman);

    kernel mean_filter = {
            {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
            3,
            9
    };

    kernel gaussian_filter = {
            {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}},
            3,
            16
    };

    kernel laplace_filter = {
            {{0, -1, 0}, {-1, 4, -1}, {0, -1, 0}},
            3,
            1
    };

    kernel high_pass_filter = {
            {{-1, -1, -1}, {-1, 9, -1}, {-1, -1, -1}},
            3,
            1
    };

    Mat low_passed_mean = filter(cameraman, mean_filter);
    Mat low_passed_gaussian = filter(cameraman, gaussian_filter);
    Mat high_passed_laplace = filter(cameraman, laplace_filter);
    Mat high_passed = filter(cameraman, high_pass_filter);

    imshow("Low pass: mean filter", low_passed_mean);
    imshow("Low pass: gaussian filter", low_passed_gaussian);
    imshow("High pass: laplace filter", high_passed_laplace);
    imshow("High pass: filter", high_passed);

    Mat frequency_domain = generic_frequency_domain_filter(cameraman);
    imshow("Frequency domain filtering", frequency_domain);

    waitKey();

    return 0;
}
