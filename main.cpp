#include <iostream>
#include <opencv2/opencv.hpp>
#include <fstream>


using namespace std;
using namespace cv;

typedef struct neighborhood_structure{
    int size;
    int* di;
    int* dj;
};

bool IsInside(Mat img, int i, int j){
    return (i >= 0 && i < img.rows)  && (j >= 0 && j < img.cols);
}

Mat dilation(Mat source, neighborhood_structure neighborhood, int no_iter){

    //TODO: Implement the dilation operation for no_iter times using the structuring element defined by
    // the neighborhood argument

    int rows = source.rows, cols = source.cols;
    Mat dst(rows, cols, CV_8UC1, Scalar(255)), aux = source;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int m = 0; m < no_iter; m++) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (aux.at<uchar>(i, j) == 0) {
                    dst.at<uchar>(i, j) = 0;

                    for (int k = 0; k < neighborhood.size; k++) {
                        int new_i = i + neighborhood.di[k];
                        int new_j = j + neighborhood.dj[k];

                        if (IsInside(aux, new_i, new_j)) {
                            dst.at<uchar>(new_i, new_j) = 0;
                        }
                    }
                }
            }
        }

        aux = dst.clone();
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return dst;

}

Mat erosion(Mat source, neighborhood_structure neighborhood, int no_iter){

    //TODO: Implement the erosion operation for no_iter times using the structuring element defined by
    // the neighborhood argument

    int rows = source.rows, cols = source.cols;
    Mat dst(rows, cols, CV_8UC1, Scalar(255)), aux = source;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int m = 0; m < no_iter; m++) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                bool background = false;

                for (int k = 0; k < neighborhood.size; k++) {
                    int new_i = i + neighborhood.di[k];
                    int new_j = j + neighborhood.dj[k];

                    if (IsInside(aux, new_i, new_j)) {
                        if (aux.at<uchar>(new_i, new_j) == 255) {
                            background = true;

                            break;
                        } else {
                            background = false;
                        }
                    }
                }

                dst.at<uchar>(i, j) = background ? 255 : 0;
            }
        }

        aux = dst;
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****


    return dst;
}

Mat opening(Mat source, neighborhood_structure neighborhood, int no_iter) {

    //TODO: Implement the opening operation for no_iter times using the structuring element defined by
    // the neighborhood argument

    Mat dst = source;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < no_iter; i++) {
        dst = dilation(erosion(dst, neighborhood, 1), neighborhood, 1);
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return dst;

}

Mat closing(Mat source, neighborhood_structure neighborhood, int no_iter) {

    //TODO: Implement the closing operation for no_iter times using the structuring element defined by
    // the neighborhood argument

    Mat dst = source;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < no_iter; i++) {
        dst = erosion(dilation(dst, neighborhood, 1), neighborhood, 1);
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return dst;
}

Mat boundary_extraction(Mat source, neighborhood_structure neighborhood) {

    //TODO: Implement the boundary extraction algorithm for no_iter times using the structuring element defined by
    // the neighborhood argument

    int rows = source.rows, cols = source.cols;
    Mat erosion_mat = erosion(source, neighborhood, 1), dst(rows, cols, CV_8UC1, Scalar(0));

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (source.at<uchar>(i, j) == erosion_mat.at<uchar>(i, j)) {
                dst.at<uchar>(i, j) = 255;
            } else {
                dst.at<uchar>(i, j) = 0;
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return dst;

}

Mat region_filling(Mat source, neighborhood_structure neighborhood) {

    //TODO: Implement the region filling algorithm for no_iter times using the structuring element defined by
    // the neighborhood argument

    Mat dst;
    int rows, cols;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return dst;

}

int main() {
    Mat source = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 7/Morphological_Op_Images/1_Dilate/balloons_bw.bmp",
                        IMREAD_GRAYSCALE);
    Mat source_boundary = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 7/Morphological_Op_Images/5_BoundaryExtraction/reg1neg1_bw.bmp",
                                 IMREAD_GRAYSCALE);

    imshow("Original Image", source);

    neighborhood_structure n8 = {
            9,
            (int[]){0,-1,-1,-1,0,1,1,1, 0},
            (int[]){1,1,0,-1,-1,-1,0,1, 0}
    };

    neighborhood_structure n4 = {
            5,
            (int[]){ 0,-1, 0, 1, 0},
            (int[]){ 1, 0,-1, 0, 0}
    };

    Mat dilation_4n = dilation(source, n4, 1);
//    imshow("Dilation 4n", dilation_4n);

    Mat dilation_8n = dilation(source, n8, 1);
//    imshow("Dilation 8n", dilation_8n);

    Mat erosion_4n = erosion(source, n4, 1);
    imshow("Erosion 4n", erosion_4n);

    Mat erosion_8n = erosion(source, n8, 1);
    imshow("Erosion 8n", erosion_8n);

    Mat opening_4n = opening(source, n4, 2);
//    imshow("Opening 4n", opening_4n);

    Mat opening_8n = opening(source, n8, 2);
//    imshow("Opening 8n", opening_8n);

    Mat closing_4n = closing(source, n4, 2);
//    imshow("Closing 4n", closing_4n);

    Mat closing_8n = closing(source, n8, 2);
//    imshow("Closing 8n", closing_8n);



    Mat boundary_extraction_mat = boundary_extraction(source_boundary, n8);
//    imshow("Boundary extraction", boundary_extraction_mat);

//    Mat region_filling_mat = region_filling(source_boundary, n8);
//    imshow("Region filling", region_filling_mat);

    waitKey();

    return 0;
}
