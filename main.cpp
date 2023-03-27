#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int n8_di[8] = {0,-1,-1, -1, 0, 1, 1, 1};
int n8_dj[8] = {1, 1, 0, -1, -1,-1, 0, 1};

int np_di[4] = { 0,-1,-1, -1};
int np_dj[4] = { -1,-1, 0, 1};

/*
 * Structure storing the label matrix in labels and
 * the number of labels in no_labels
 */
struct labels{
    Mat labels;
    int no_labels;
};

Mat color_labels(labels labels_str){

    /*
     * This method will generate a number of no_labels colors and
     * generate a color image containing each connected component displayed in a different color
     */

    Mat labels = labels_str.labels;
    int rows = labels.rows, cols = labels.cols, no_labels = labels_str.no_labels;
    Mat result(rows, cols, CV_8UC3, Vec3b(255, 255, 255));

    Vec3b* colors = static_cast<Vec3b *>(calloc(255, sizeof(Vec3b)));

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uchar label = labels.at<uchar>(i, j);
            if (label != 0) {
                if (colors[label][0] == 0 && colors[label][1] == 0 && colors[label][2] == 0) {
                    Vec3b color(rand() % 255, rand() % 255, rand() % 255);
                    colors[label] = color;
                }

                result.at<Vec3b>(i, j) = colors[label];
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return result;

}

labels BFS_labeling(Mat source){

    /*
     * This method will implement the BFS labeling algorithm
     * Hint:
     *  Use the Point structure(or a similar one) to store the coordinates in a queue
     *  You can use queue from C++ with its specific operations (push, pop, empty, front)
     */
    int rows = source.rows, cols = source.cols;
    Mat labels(rows, cols, CV_8UC1, Scalar(0));
    int label = 0;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (source.at<uchar>(i, j) == 0 && labels.at<uchar>(i, j) == 0) {
                label++;
                labels.at<uchar>(i, j) = label;

                queue<Point> q;
                q.push(Point(i, j));

                while (!q.empty()) {
                    Point head = q.front();
                    q.pop();

                    for (int k = 0; k < 8; k++) {
                        for (int l = 0; l < 8; l++) {
                            int x = head.x + n8_di[k];
                            int y = head.y + n8_dj[l];

                            if (source.at<uchar>(x, y) == 0 && labels.at<uchar>(x,y) == 0) {
                                labels.at<uchar>(x, y) = label;
                                q.push(Point(x, y));
                            }
                        }
                    }
                }
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return {labels, label};
}

labels Two_pass_labeling(Mat source){
    Mat labels;
    int rows, cols;
    int newlabel, label;
    /*
     * This method will implement the two pass labeling algorithm
     * Hint:
     *  Use the vector structure from C++(actually you need a vector of vectors and a simple one, check out the lab works)
     *  You can use queue from C++ with its specific operations (push, pop, empty, front)
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****



    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****


    return {labels, newlabel};
}
int main() {
    Mat source = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 5/PI-L5/letters.bmp",
                        IMREAD_GRAYSCALE);

    imshow("Original Image", source);

    labels bfsLabels = BFS_labeling(source);
    Mat result_bfs = color_labels(bfsLabels);
    imshow("BFS", result_bfs);

//    labels two_pass_label = Two_pass_labeling(source);
//    Mat result_two_pass = color_labels(two_pass_label);
//    imshow("Two pass", result_two_pass);

    waitKey();

    return 0;
}
