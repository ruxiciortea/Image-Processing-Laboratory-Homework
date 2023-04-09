#include <opencv2/opencv.hpp>
#include <fstream>
#include <iostream>

using namespace std;
using namespace cv;

int n8_di[8] = {0,-1,-1, -1, 0, 1, 1, 1};
int n8_dj[8] = {1, 1, 0, -1, -1,-1, 0, 1};

Point find_P_0(Mat source){
    /*
     * Find the initial point of the contour and return it
     */
    Point P_0;
    int rows = source.rows, cols = source.cols;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (source.at<uchar>(i, j) == 0) {
                P_0 = Point(i, j);
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return P_0;
}

typedef struct contour{
    vector<Point> border;
    vector<int> dir_vector;
};

contour extract_contour(Mat source, Point P_0){

    /*
     * Use the border tracing algorithm in order to extract the contour
     * Save it as a vector of points and a vector of directions
     */

    int dir = 7;
    Point P_current = P_0;

    vector<Point> border(source.rows * source.cols, Point(0,0));
    vector<int> dir_vector(source.rows * source.cols, 0);

    border[0] = P_current;
    dir_vector[0] = dir;

    int border_size = 1, dir_size = 1;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    while (border_size <= 2 || (border[0] != border[border_size - 2] && border[1] != border[border_size - 1]) || 0) {
        if (dir % 2 == 0) {
            dir = (dir + 7) % 8;
        } else {
            dir = (dir + 6) % 8;
        }

        for (int i = 0; i < 8 ; i++) {
            int current_dir = dir + i;

            int new_x = P_current.x + n8_dj[current_dir % 8];
            int new_y = P_current.y + n8_di[current_dir % 8];

            uchar neighbour = source.at<uchar>(new_y, new_x);

            if (neighbour == 0) {
                P_current = Point(new_x, new_y);
                border[border_size++] = P_current;

                dir = current_dir % 8;
                dir_vector[dir_size++] = dir;

                break;
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return {border, dir_vector};
}

Mat draw_contour(contour cnt, Mat source){

    /*
     * Draw the contour using the border variable from cnt structure
     */

    int rows = source.rows, cols = source.cols;
    Mat dst(rows, cols, CV_8UC1, Scalar(255));

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < cnt.border.size(); i++) {
        Point current_point = cnt.border[i];

        dst.at<uchar>(current_point.y, current_point.x) = 0;
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****


    return dst;
}

void print_AC_DC_chain_codes(contour cnt){
    /*
     * Print the AC and DC chain codes
     * Hint: You have to compute the DC one
     */

    printf("The AC vector is: ");

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    printf("\nThe DC vector is: ");

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

Mat contour_reconstruction(FILE *pf, Mat background) {

    /*
     * From the file read the chain code and reconstruct the image
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return background;

}


int main() {
    Mat source = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 6/PI-L6/star_R90.bmp",
                        IMREAD_GRAYSCALE);

    imshow("Original Image", source);

    Point P_0 = find_P_0(source);
    contour cnt = extract_contour(source, P_0);
    Mat mat_cnt = draw_contour(cnt, source);

    imshow("Contour", mat_cnt);

//    print_AC_DC_chain_codes(cnt);
//
//    FILE *pf;
//    pf = fopen("/Users/ruxiciortea/Desktop/IP/Labs/Lab 6/PI-L6/reconstruct.txt", "r");
//    Mat background = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 6/PI-L6/gray_background.bmp", IMREAD_GRAYSCALE);
//
//    Mat reconstruction = contour_reconstruction(pf, background);
//
//    imshow("Reconstruction", reconstruction);

    waitKey();

    return 0;
}
