//
// Created by GuoweiWu on 8/21/18.
//

#include <iostream>
#include <thread>
#include <chrono>
#include "commonUtils.h"
#include "MatrixTranspose.h"
#include "MatrixMulplication.h"

int main() {
    transposeMatrix();     // compare three different transpose methods

    matrixMulplication();  // compare three different multiply methods
}