//
// Created by GuoweiWu on 8/21/18.
//
#include <iostream>
#include <chrono>

#define M 1600
#define N 1600
#define K 1600
#define TILE 16
#define WIDTH 3

void print(int* matrix, int row, int col);
void bruteForceMultiply(int* A, int* B, int* C);
void tiledMultiply(int* A, int* B, int* C, int tileSize);
void compareArrays(int* m1, int *m2, int row, int col);

int main() {

    // initialize random generator
    srand(2018);

    // declare matrix
    int* A = (int*)malloc(sizeof(int) * M * N);
    int* B = (int*)malloc(sizeof(int) * N * K);
    int* C_BruteForce = (int*)malloc(sizeof(int) * M * K);
    int* C_Tiled = (int*)malloc(sizeof(int) * M * K);
    memset(C_Tiled, 0, M*K*sizeof(int));

    // Generate Matrix using random number in range [0,99]
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i * N + j] = rand() % 10;
        }
    }
//    print(A, M, N);

    // Generate Matrix using random number in range [0,99]
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < K; ++j) {
            B[i * K + j] = rand() % 10;
        }
    }
//    print(B, N, K);

    auto start_bf = std::chrono::high_resolution_clock::now();
    bruteForceMultiply(A, B, C_BruteForce);
    auto end_bf = std::chrono::high_resolution_clock::now();
    auto duration_bf = std::chrono::duration_cast<std::chrono::microseconds>(end_bf - start_bf);
//    print(C_BruteForce, M, K);

    auto start_tile = std::chrono::high_resolution_clock::now();
    tiledMultiply(A, B, C_Tiled, TILE);
    auto end_tile = std::chrono::high_resolution_clock::now();
    auto duration_tile = std::chrono::duration_cast<std::chrono::microseconds>(end_tile - start_tile);
//    print(C_Tiled, M, K);

    compareArrays(C_BruteForce, C_Tiled, M, K);

    printf("exec time for brute force: %d milliseconds\n", duration_bf.count()/1000);
    printf("exec time for tile method: %d milliseconds\n", duration_tile.count()/1000);
}

// utility to print a matrix
void print(int* matrix, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%*d ", WIDTH, matrix[i * col + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void bruteForceMultiply(int *A, int *B, int *C) {
    for (int row = 0; row < M; row++) {
        for (int col = 0; col < K; col++) {
            int cur = 0;
            for (int inner = 0; inner < N; inner++) {
                cur += A[row * N + inner] * B[inner * K + col];
            }
            C[row * K + col] = cur;
        }
    }
}

void tiledMultiply(int* A, int* B, int* C, int tileLen) {
    for (int i = 0; i + tileLen - 1 < M; i += tileLen) {
        for (int j = 0; j + tileLen - 1 < N; j += tileLen) {
            for (int h = 0; h + tileLen - 1 < K; h += tileLen) {
                // calculating a square tile unit
                for (int r = i; r < i + tileLen; r++) {
                    for (int inner = j; inner < j + tileLen; inner++) {
                        for (int c = h; c < h + tileLen; c++) {
                            C[r * K + c] += A[r * N + inner] * B[inner * K + c];
                        }
                    }
                }
            }
        }
    }
}

void compareArrays(int* m1, int *m2, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (m1[i * col + j] != m2[i * col + j]) printf(" ====== Not equal======\n");
        }
    }
    printf(" ====== equal======\n");
}