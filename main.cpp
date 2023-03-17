#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

#define HISTOGRAM_SIZE 256

bool IsInside(Mat img, int i, int j){
    /*
    * Implement a function called isInside(img, i, j) which checks if the position indicated by
    * the pair (i,j) (row, column) is inside the image img.
    */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    return (i >= 0 && i <= img.rows)  && (j >= 0 && j <= img.cols);

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

typedef struct grayscale_mapping{
    uchar* grayscale_values; //hold the grayscale values after thresholding
    uchar count_grayscale_values; //hold the number grayscale values after thresholding
};

int* compute_histogram(Mat source, int histogram_bins){

    /*
    * Compute  the  histogram  for  a  given  grayscale  image (in  an  array  of  integers  having dimension 256)
    */

    int rows = source.rows;
    int cols = source.cols;
    int* histogram = new int[HISTOGRAM_SIZE];

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < HISTOGRAM_SIZE; ++i) {
        histogram[i] = 0;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            unsigned char pixel = source.at<uchar>(i, j);
            histogram[pixel]++;
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return histogram;

}

int* compute_histogram_custom(Mat source, int histogram_bins){

    /*
     * Compute the histogram for a given number of bins m≤ 256.
     */

    int rows = source.rows;
    int cols = source.cols;
    int* histogram = new int[HISTOGRAM_SIZE];
    int bins_in_section = HISTOGRAM_SIZE / histogram_bins;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < HISTOGRAM_SIZE; ++i) {
        histogram[i] = 0;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            unsigned char pixel = source.at<uchar>(i, j);
            int index = pixel / bins_in_section;
            histogram[index]++;
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return histogram;

}

float* compute_pdf(int* histogram, Mat source){
    /*
     *Compute the PDF (in an array of floats of dimension 256)
     */

    int rows = source.rows;
    int cols = source.cols;
    int img_size = rows * cols;
    int no_grayscale_values = HISTOGRAM_SIZE;
    float* pdf = new float[HISTOGRAM_SIZE];

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < HISTOGRAM_SIZE; ++i) {
        pdf[i] = 0;
    }

    for (int i = 0; i < no_grayscale_values; i++) {
        pdf[i] = (float)histogram[i] / img_size;
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return pdf;

}

void showHistogram(const string& name, int* hist, const int  hist_cols, const int hist_height){
    /*
     * Hint: Look in the lab work
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

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

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

grayscale_mapping multi_level_thresholding(Mat source, int wh, float th, float* pdf){
    /*
     * Implement the multilevel thresholding algorithm from section 3.3.
     */

    grayscale_mapping map;
    int window_width = 2 * wh + 1;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    map.count_grayscale_values = 1;
    map.grayscale_values = new uchar[HISTOGRAM_SIZE];

    for (int i = 0; i < HISTOGRAM_SIZE; ++i) {
        map.grayscale_values[i] = 0;
    }

    for (int k = wh; k <= 255 - wh; k++) {
        float average = 0;
        bool is_max = true;

        for (int i = k - wh; i <= k + wh; i++) {
            average += pdf[i];

            if (pdf[k] < pdf[i]) {
                is_max = false;
            }
        }

        average /= window_width;

        if (pdf[k] > average + th) {
            if (is_max) {
                map.grayscale_values[map.count_grayscale_values++] = k;
            }
        }
    }

    map.grayscale_values[map.count_grayscale_values] = 255;

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return map;

}

uchar find_closest_histogram_maximum(uchar old_pixel, grayscale_mapping gray_map){

    /*
     * Find the corresponding quantized value to map a pixel
     * Hint: Look in the gray_map and find out the value that resides at index argmin of the distance between old_pixel
     *      and the values in gray_map
     */

    uchar new_grayscale_value;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    int min_dist = 255;

    for (int i = 0; i < gray_map.count_grayscale_values; i++) {
        int distance = abs(gray_map.grayscale_values[i] - old_pixel);

        if (distance < min_dist) {
            min_dist = distance;
            new_grayscale_value = gray_map.grayscale_values[i];
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****


    return new_grayscale_value;
}

Mat draw_multi_thresholding(Mat source, grayscale_mapping grayscale_map){

    /*
     * Draw the new multi level threshold image by mapping each pixel to the corresponding quantized values
     * Hint: Look in the grayscale_map structure for all the obtained grayscale values and for each pixel in the
     *      source image assign the correct value. You may use the find_closest_histogram_maximum function
     */

    int rows = source.rows;
    int cols = source.cols;
    Mat result(rows, cols, CV_8UC1);

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.at<uchar>(i, j) = find_closest_histogram_maximum(source.at<uchar>(i, j), grayscale_map);
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****


    return result;
}

uchar update_pixel_floyd_steinberg_dithering(uchar pixel_value, int value){
    /*
     * Update the value of a pixel in the floyd_steinberg alg.
     * Take care of the values bellow 0 or above 255. Clamp them.
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    int rez = pixel_value + value;

    if (rez < 0) {
        return 0;
    }

    if (rez > 255) {
        return 255;
    }

    return rez;

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

Mat floyd_steinberg_dithering(Mat source, grayscale_mapping grayscale_map){

    /*
     * Enhance  the  multilevel  thresholding  algorithm  using  the  Floyd-Steinberg  dithering from section 3.4.
     * Hint: Use the update_pixel_floyd_steinberg_dithering when spreading the error
     */

    int rows = source.rows;
    int cols = source.cols;
    Mat result(rows, cols, CV_8UC1);

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows -1; i++) {
        for (int j = 0; j < cols - 1; j++) {
            uchar old_pixel = source.at<uchar>(i,j);
            uchar new_pixel = find_closest_histogram_maximum(old_pixel, grayscale_map);
            result.at<uchar>(i, j) = new_pixel;

            uchar error = abs(old_pixel - new_pixel);

            if (IsInside(source, i, j + 1)) {
                uchar pixel_value = source.at<uchar>(i, j + 1);
                result.at<uchar>(i, j + 1) = update_pixel_floyd_steinberg_dithering(pixel_value, 7 * error / 16);
            }

            if (IsInside(source, i + 1, j - 1)) {
                uchar pixel_value = source.at<uchar>(i + 1, j - 1);
                result.at<uchar>(i + 1, j - 1) = update_pixel_floyd_steinberg_dithering(pixel_value, 3 * error / 16);
            }

            if (IsInside(source, i + 1, j)) {
                uchar pixel_value = source.at<uchar>(i + 1, j);
                result.at<uchar>(i + 1, j) = update_pixel_floyd_steinberg_dithering(pixel_value, 5 * error / 16);
            }

            if (IsInside(source, i + 1, j + 1)) {
                uchar pixel_value = source.at<uchar>(i + 1, j + 1);
                result.at<uchar>(i + 1, j + 1) = update_pixel_floyd_steinberg_dithering(pixel_value, error / 16);
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

    return result;
}

int main() {
    Mat cameraman = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 3/PI-L3/cameraman.bmp",
                           IMREAD_GRAYSCALE);
    Mat saturn = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 3/PI-L3/saturn.bmp",
                        IMREAD_GRAYSCALE);

    imshow("Cameraman original", cameraman);
    imshow("Saturn original", saturn);

    int* histogram_cameraman = compute_histogram(cameraman, 256);
    float* pdf_cameraman = compute_pdf(histogram_cameraman, cameraman);

    int* histogram_saturn = compute_histogram(saturn, 256);
    float* pdf_saturn = compute_pdf(histogram_saturn, saturn);

    printf("Some histogram values are: ");
    for(int i = 50; i < 56; i++){
        printf("%d ", histogram_cameraman[i]);
    }
    printf("\n");

    printf("Some pdf values are: ");
    for(int i = 50; i < 56; i++){
        printf("%f ", pdf_cameraman[i]);
    }

    showHistogram("Histogram", histogram_cameraman, 256, 100);

    int* histogram_custom = compute_histogram_custom(cameraman, 40);
    showHistogram("Histogram reduced bins", histogram_custom, 40, 100);

    delete [] histogram_cameraman;
    delete [] histogram_saturn;
    delete [] histogram_custom;

    grayscale_mapping grayscale_map_saturn = multi_level_thresholding(saturn, 5, 0.0003, pdf_saturn);
    grayscale_mapping grayscale_map_cameraman = multi_level_thresholding(cameraman, 5, 0.0003, pdf_cameraman);

    delete [] pdf_cameraman;
    delete [] pdf_saturn;

    Mat image_multi_threshold_cameraman = draw_multi_thresholding(cameraman, grayscale_map_cameraman);
    imshow("Multi level threshold cameraman", image_multi_threshold_cameraman);

    Mat fsd_cameraman = floyd_steinberg_dithering(cameraman, grayscale_map_cameraman);
    imshow("Floyd Steinberg Dithering cameraman", fsd_cameraman);

    Mat image_multi_threshold_saturn = draw_multi_thresholding(saturn, grayscale_map_saturn);
    imshow("Multi level threshold saturn", image_multi_threshold_saturn);

    Mat fsd_saturn = floyd_steinberg_dithering(saturn, grayscale_map_saturn);
    imshow("Floyd Steinberg Dithering saturn", fsd_saturn);

    waitKey(0);
    return 0;
}
