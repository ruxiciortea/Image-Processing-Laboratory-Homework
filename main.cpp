#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

typedef struct image_channels_bgr{
    Mat B;
    Mat G;
    Mat R;
};

typedef struct image_channels_hsv{
    Mat H;
    Mat S;
    Mat V;
};

image_channels_bgr break_channels(Mat source) {

/*
* Create a function that will copy the R, G and B channels of a color,
* RGB image (CV_8UC3 type) into three matrices of type CV_8UC1 (grayscale images).
* Return the three matrices in an image_channels_bgr structure
*
* Inputs:
*  source: the source image(BGR)
* Variables:
*  rows: number of rows of the source matrix
*  cols: number of cols of the source matrix
*  B, G, R: Matrices that will store each a color channel
*  bgr_channels: structure of type image_channels_bgr that will return three channels
*/

    int rows = source.rows, cols = source.cols;
    Mat B(rows, cols, CV_8UC1), G(rows, cols, CV_8UC1), R(rows, cols, CV_8UC1);
    image_channels_bgr bgr_channels;

//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Vec3b pixel = source.at< Vec3b>(i,j);
            B.at<uchar>(i, j) = pixel[0];
            G.at<uchar>(i, j) = pixel[1];
            R.at<uchar>(i, j) = pixel[2];
        }
    }

    bgr_channels.B = B;
    bgr_channels.G = G;
    bgr_channels.R = R;

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return bgr_channels;

}

void display_channels(image_channels_bgr bgr_channels){
/*
* Display each channel in a different window
* Do not put here the waitKey() try to use only one at the end of the lab
* Input:
*  bgr_channels: structure of type image_channels_bgr that will store B, G, R channels
*
*/
//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    imshow("Blue channel", bgr_channels.B);
    imshow("Green channel", bgr_channels.G);
    imshow("Red channel", bgr_channels.R);

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

Mat bgr_2_grayscale(Mat source){
/*
* Create a function that will convert a color RGB image (CV_8UC3 type) to a grayscale image (CV_8UC1),
* and return the result image
* Inputs: source: the source matrix(BGR)
* Variables:
*  rows: number of rows of the source matrix
*  cols: number of cols of the source matrix
*  grayscale_image: The grayscale image that you will obtain and return
*/
    int rows = source.rows, cols = source.cols;
    Mat grayscale_image(rows, cols, CV_8UC1);

//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Vec3b pixel = source.at< Vec3b>(i,j);
            grayscale_image.at<uchar>(i, j) = (pixel[0] + pixel[1] + pixel[2]) / 3;
        }
    }

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return grayscale_image;

}

Mat grayscale_2_binary(Mat source, int threshold){
/*
* Create a function for converting from grayscale to black and white (binary), using (2.2).
* Test the operation on multiple images, and using multiple thresholds.
* Inputs:
*    source: grayscale image
*    threshold: the threshold you are going to use to perform the binarization
*  Variables:
*    rows: number of rows of the source matrix
*    cols: number of cols of the source matrix
*    binary: the resulted binarized image
*/

    int rows = source.rows, cols = source.cols;
    Mat binary(rows, cols, CV_8UC1);;

//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

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

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return binary;
}

image_channels_hsv bgr_2_hsv(image_channels_bgr bgr_channels){
/*
* Create a function that will compute the H, S and V values from the R, G, B channels of
* an image, using the equations from 2.6. Store each value (H, S, V) in a CV_8UC1 matrix that will be
* stored in an image_channels_hsv struct.
* Inputs:
*    bgr_channels: structure that stores the B, G, R channels of an image
* Variables:
*    rows: number of rows of the source matrix
*    cols: number of cols of the source matrix
*    H, S, V: matrices that will store the values of the 3 different channels (Pay attention to the type of elements
*       that are stored in these matrices.
*    hsv_channels: structure that will store the H, S, V channels
*/

    int rows = bgr_channels.R.rows, cols = bgr_channels.R.cols;
    Mat H(rows, cols, CV_32FC1), S(rows, cols, CV_32FC1), V(rows, cols, CV_32FC1);
    Mat r(rows, cols, CV_32FC1), g(rows, cols, CV_32FC1), b(rows, cols, CV_32FC1);
    image_channels_hsv hsv_channels;

//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float r_local = (float)bgr_channels.R.at<uchar>(i, j) / 255;
            float g_local = (float)bgr_channels.G.at<uchar>(i, j) / 255;
            float b_local = (float)bgr_channels.B.at<uchar>(i, j) / 255;

            r.at<float>(i, j) = r_local;
            g.at<float>(i, j) = g_local;
            b.at<float>(i, j) = b_local;

            float M = max(max(r_local, g_local), b_local);
            float m = min(min(r_local, g_local), b_local);

            float C = M - m;

            // value
            V.at<float>(i, j) = M;

            // saturation
            if (V.at<float>(i, j) != 0) {
                S.at<float>(i, j) = C / V.at<float>(i, j);
            } else {
                S.at<float>(i, j) = 0;
            }

            // hue
            if (C != 0) {
                if (M == r_local) {
                    H.at<float>(i, j) = 60 * (g_local - b_local) / C;
                } else if (M == g_local) {
                    H.at<float>(i, j) = 120 + 60 * (b_local - r_local) / C;
                } else if (M == b_local) {
                    H.at<float>(i, j) = 240 + 60 * (r_local - g_local) / C;
                } else { // greyscale
                    H.at<float>(i, j) = 0;
                }

                if (H.at<float>(i, j) < 0) {
                    H.at<float>(i, j) = H.at<float>(i, j) + 360;
                }
            }
        }
    }

    hsv_channels.H = H;
    hsv_channels.S = S;
    hsv_channels.V = V;

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return hsv_channels;
}

void display_hsv_channels(image_channels_hsv hsv_channels){

/*
* Display the three channels
* Inputs:
*    hsv_channels: structure that stores the H, S, V channels of an image.
*      In order to display them don't forget to normalize them accordingly
*  Variables:
*    rows: number of rows of the source matrix
*    cols: number of cols of the source matrix
*    H_norm, S_norm, V_norm: Normalized matrices.
*/

    int rows = hsv_channels.H.rows, cols = hsv_channels.H.cols;
    Mat H_norm(rows, cols, CV_8UC1), S_norm(rows, cols, CV_8UC1), V_norm(rows, cols, CV_8UC1);

//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float h = hsv_channels.H.at<float>(i, j);
            float s = hsv_channels.S.at<float>(i, j);
            float v = hsv_channels.V.at<float>(i, j);

            H_norm.at<uchar>(i, j) = h * 255 / 360;
            S_norm.at<uchar>(i, j) = s * 255;
            V_norm.at<uchar>(i, j) = v * 255;
        }
    }

    imshow("H channel", H_norm);
    imshow("S channel", S_norm);
    imshow("V channel", V_norm);

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

bool IsInside(Mat img, int i, int j){
/*
* Implement a function called isInside(img, i, j) which checks if the position indicated by
* the pair (i,j) (row, column) is inside the image img.
*/

//*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    return (i >= 0 && i <= img.rows)  && (j >= 0 && j <= img.cols);

//*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}


int main() {
    Mat image = imread("/Users/ruxiciortea/Desktop/IP/Lab 2/PI-L2/flowers_24bits.bmp",
                       IMREAD_COLOR);

    imshow("Original image", image);

    image_channels_bgr bgr_channels = break_channels(image);
    display_channels(bgr_channels);

    Mat greyscale_image = bgr_2_grayscale(image);
    imshow("Grayscale image", greyscale_image);

//Try with 30, 60, 180, 220
    Mat binary_image = grayscale_2_binary(greyscale_image, 60);
    imshow("Binary image", binary_image);

    image_channels_hsv hsv_channels = bgr_2_hsv(bgr_channels);
    display_hsv_channels(hsv_channels);

    bool is_inside = IsInside(greyscale_image, 300, 50);
    cout << "Rows in image: " << greyscale_image.rows << endl;
    cout << "Columns in image: " << greyscale_image.cols << endl;
    cout << "Is the given pixel inside the image? " << is_inside;

    waitKey(0);
    return 0;
}
