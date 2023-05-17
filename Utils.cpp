//
// Created by Ruxandra Ciortea on 17.05.2023.
//

#include "Utils.h"

float** initMatrix(int rows, int cols) {
    float **matrix = (float**)calloc(rows, sizeof(float*));

    for (int i = 0; i < rows; i++) {
        matrix[i] = (float*)calloc(cols, sizeof(float));
    }

    return matrix;
}

void printMatrix(int rows, int cols, float **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << matrix[i][j] << " ";
        }

        cout << "\n";
    }

    cout << "\n";
}

int findRegion(double angle) {
    if ((angle >= 67.5 && angle < 112.5) || (angle >= 247.5 && angle < 292.5)) {
        return 0;
    }

    if ((angle >= 22.5 && angle < 67.5) || (angle >= 202.5 && angle < 247.5)) {
        return 1;
    }

    if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle < 202.5) || (angle >= 337.5 && angle <= 360)) {
        return 2;
    }

    if ((angle >= 112.5 && angle < 157.5) || (angle >= 292.5 && angle < 337.5)) {
        return 3;
    }

    return -1;
}

void findPointsBasedOnRegion(int region, Point &p1, Point &p2) {
    switch (region) {
        case 0:
            p1 = Point(0, -1);
            p2 = Point(0, 1);

            break;
        case 1:
            p1 = Point(1, -1);
            p2 = Point(-1, 1);

            break;
        case 2:
            p1 = Point(1, 0);
            p2 = Point(-1, 0);

            break;
        case 3:
            p1 = Point(-1, -1);
            p2 = Point(1, 1);

            break;
    }
}

Kernel initKernel(vector<int> values, int size) {
    float **matrix = initMatrix(size, size);
    int index = 0;
    float sum = 0;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i][j] = values[index];
            sum += values[index];
            index++;
        }
    }

    float meanValue = sum == 0 ? 1 : sum;

    return {matrix, size, meanValue};
}
