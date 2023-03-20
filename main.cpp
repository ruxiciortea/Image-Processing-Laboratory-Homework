#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

#define PI 3.14

typedef struct perimeter{
    Mat contour;
    int length;
};

typedef struct circumscribed_rectangle_coord{
    int c_min;
    int c_max;
    int r_min;
    int r_max;
};

typedef struct elongation_axis{ //openCV points are X,Y (column, rows)
    Point p1;
    Point p2;
};

int y_coord_neighborhood[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
int x_coord_neighborhood[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

bool compare_pixels(Vec3b pixel_1, Vec3b pixel_2){
    //This method will compare if two pixels are equal or not(if it's the same color or not)
    //Return true if pixel_1 == pixel_2

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    return pixel_1[0] == pixel_2[0] && pixel_1[1] == pixel_2[1] && pixel_1[2] == pixel_2[2];
    
    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) *****

}

Mat get_object_instance(Mat source, Vec3b color){
    /*
     * This method will save the selected object, in a different matrix, in a binary format(0 and 255)
     */

    int rows = source.rows;
    int cols = source.cols;

    Mat result(rows, cols, CV_8UC1, Scalar(255));
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (compare_pixels(source.at<Vec3b>(i, j), color)) {
                result.at<uchar>(i, j) = 0;
            }
        }
    }
    
    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return result;
}

perimeter naive_perimeter(Mat binary_object){

    /*
     * This method will compute the perimeter and save the contour in a perimeter structure
     * that will store the two components
     */
    int rows = binary_object.rows;
    int cols = binary_object.cols;

    perimeter object_perimeter;
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    object_perimeter.contour = Mat(rows, cols, CV_8UC1, Scalar(255));
    object_perimeter.length = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (binary_object.at<uchar>(i, j) == 0) {
                bool found = false;

                for (int k = 0; k < 8; k++) {
                    for (int l = 0; l < 8; l++) {
                        if (binary_object.at<uchar>(i + x_coord_neighborhood[k], j + y_coord_neighborhood[l]) == 255
                                && !found) {
                            object_perimeter.contour.at<uchar>(i, j) = 0;
                            object_perimeter.length++;

                            found = true;
                        }
                    }
                }
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****


    return object_perimeter;

}

int compute_area(Mat binary_object){
    /*
     * This method will compute the object area and return it
     */
    int rows = binary_object.rows;
    int cols = binary_object.cols;

    int area = 0;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (binary_object.at<uchar>(i, j) == 0) {
                area++;
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return area;
}

Point compute_center_of_mass(Mat binary_object){
    /*
     * This method will compute the center of the mass and return it in a Point structure
     * Hint: Pay attention at the x and y notation
     */
    int rows = binary_object.rows;
    int cols = binary_object.cols;

    Point center_mass;
    int rSum = 0;
    int cSum = 0;
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (binary_object.at<uchar>(i, j) == 0) {
                rSum += i;
                cSum += j;
            }
        }
    }

    int area = compute_area(binary_object);

    center_mass.x = cSum / area;
    center_mass.y = rSum / area;

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return center_mass;

}

Mat display_center_of_mass(Point center_of_mass, Mat source){

    /*
     * This method will display on the source image the center_of_mass
     * Hint: Use the circle method from OpenCv, clone the source
     */
    Mat result = source.clone();
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    circle(result, center_of_mass, 3, Scalar(255));

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return result;

}

circumscribed_rectangle_coord compute_circumscribed_rectangle_coord(Mat binary_object){
    /*
     * This method will compute the points that form the circumscribed rectangle
     * Hint:
     * c -> column
     * r -> row
     */

    int rows = binary_object.rows;
    int cols = binary_object.cols;

    circumscribed_rectangle_coord coords;
    int c_min = cols, c_max = 0, r_min = rows, r_max = 0;
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (binary_object.at<uchar>(i, j) == 0) {
                if (j < c_min) {
                    c_min = j;
                }

                if (j < c_max) {
                    c_max = j;
                }

                if (i < r_min) {
                    r_min = i;
                }

                if (i < r_max) {
                    r_max = i;
                }
            }
        }
    }

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    coords = {c_min, c_max, r_min, r_max};

    return coords;
}

float compute_aspect_ratio(circumscribed_rectangle_coord coord){
    /*
     * This method will compute the aspect ratio and will return it
     */
    float R;
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    R = coord.c_max - coord.c_min + 1 / coord.r_max - coord.r_min + 1;

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return R;
}

float compute_thinness_ratio(int area, int perimeter){
    /*
     * This method will compute the thinness ratio and will return it
     */
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
    float perimeterSquared = perimeter * perimeter;
    float fraction = area / perimeterSquared;
    float T = 4 * 3.14 * (fraction);

    return T;

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

}

float compute_axis_of_elongation_angle(Point center_of_mass, Mat binary_object){
    /*
     * This method will compute angle corresponding to the axis of elongation of a labeled object
     */
    int rows = binary_object.rows;
    int cols = binary_object.cols;

    float phi;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****


    return phi;


}

elongation_axis compute_elongation_axis_points(float angle, Point center_of_mass, circumscribed_rectangle_coord coord){
    /*
     * This method will compute the two point that are required to draw the axis of elongation and return them using the
     * elongation_axis structure
     */

    Point A, B;

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return {A, B};
}

void draw_elongation_axis(Mat source, elongation_axis axis){
    /*
     * This method will draw and display the elongation axis.
     * Hint: Use line from OpenCv, clone the source
     */

    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****


}

Mat horizontal_projection(Mat binary_image, circumscribed_rectangle_coord coord){
    /*
     * This method will compute the matrix representing the horizontal projection and return it
     */
    Mat horizontal_projection;
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return horizontal_projection;
}

Mat vertical_projection(Mat binary_image, circumscribed_rectangle_coord coord){
    /*
     * This method will compute the matrix representing the vertical projection and return it
     */
    Mat vertical_projection;
    //*****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

    //*****END OF YOUR CODE(DO NOT DELETE / MODIFY THIS LINE) ****

    return vertical_projection;
}


void geom_features(int event, int x, int y, int flags, void* param){
    Mat source = *(Mat*)param;
    Vec3b color_selected = {0, 0, 0};
    Mat binary_object;
    perimeter object_perimeter;
    Point center_of_mass;
    Mat center_of_mass_image;
    circumscribed_rectangle_coord circumscribed_coord;
    float thinness_ratio;
    float aspect_ratio;
    int area;
    float phi;
    elongation_axis axis_points;


    if (event == EVENT_LBUTTONDOWN){
        color_selected = source.at<Vec3b>(y, x);
        binary_object = get_object_instance(source, color_selected);
        imshow("Binary Object", binary_object);

        object_perimeter = naive_perimeter(binary_object);
        imshow("Contour", object_perimeter.contour);
        printf("The perimeter has length %d\n", object_perimeter.length);

        area = compute_area(binary_object);
        printf("The area is %d\n", area);

        center_of_mass = compute_center_of_mass(binary_object);
        center_of_mass_image = display_center_of_mass(center_of_mass, source);
        imshow("Center of mass", center_of_mass_image);

//        circumscribed_coord = compute_circumscribed_rectangle_coord(binary_object);
//        thinness_ratio = compute_aspect_ratio(circumscribed_coord);
//        printf("The aspect ratio is %.2f\n", thinness_ratio);

//        aspect_ratio = compute_thinness_ratio(area, object_perimeter.length);
//        printf("The thinness ratio is %.2f \n", aspect_ratio);

        aspect_ratio = compute_aspect_ratio(circumscribed_coord);
        printf("The aspect ratio is %.2f \n", aspect_ratio);

        circumscribed_coord = compute_circumscribed_rectangle_coord(binary_object);
        thinness_ratio = compute_thinness_ratio(area, object_perimeter.length);
        printf("The thinness ratio is %.2f\n", thinness_ratio);

//        phi = compute_axis_of_elongation_angle(center_of_mass, binary_object);
//        printf("The angle phi is %.2f", phi);
//        axis_points = compute_elongation_axis_points(phi, center_of_mass, circumscribed_coord);
//        draw_elongation_axis(source, axis_points);
    }
}

int main() {
    Mat image = imread("/Users/ruxiciortea/Desktop/IP/Labs/Lab 4/PI-L4/trasaturi_geom.bmp",
                       IMREAD_COLOR);

    namedWindow("Lab4", 1);

    //set the callback function for any mouse event
    setMouseCallback("Lab4", geom_features, &image);

    //show the image
    imshow("Lab4", image);

    waitKey(0);

    return 0;
}
