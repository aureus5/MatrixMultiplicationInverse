//
// Created by GuoweiWu on 8/21/18.
//
#include <iostream>
#include <thread>
#include <chrono>
#include "commonUtils.h"
#include "MatrixTranspose.h"

#define THREAD_NUM 8
//#define WIDTH 3

//void print(int* matrix, int row, int col);
void generateMatrixUsingRandomNumbers(int* matrix, int dim1, int dim2, int range);
long long getExecTime(int* A, int* B, int* C, int M, int N, int K, int TILE, int algorithm);
void bruteForceMultiply(int* A, int* B, int* C, int M, int N, int K);
void tiledMultiply(int* A, int* B, int* C, int M, int N, int K, int TILE);
void tiledMultiply_multiThread(int* A, int* B, int* C, int M, int N, int K, int TILE);
void residualBlockThread(int* A, int* B, int* C, int M, int N, int K, int TILE);
void tilingPartThreading(int* A, int* B, int* C, int M, int N, int K, int TILE, int thread_id);
void compareArrays(int* m1, int* m2, int* m3, int row, int col);

int main() {
    transposeMatrix(3425, 5234, 20);

//
//
//    // initialize random generator
//    srand(2018);
//    int range = 100; // matrix elements are randoms numbers in the the range[0, range-1], inclusive
//
//    int parameters[][4] = {{45, 35, 29, 20}, {135, 175, 124, 20},{374, 567, 421, 20},{635, 783, 323, 20},
//                           {834, 895,798, 20}, {1376, 1487, 1533, 20}, {1687, 1589, 1863, 20}};
//
//    for (int i = 0; i < sizeof(parameters) / sizeof(parameters[0]); i++) {
//        int M = parameters[i][0];
//        int N = parameters[i][1];
//        int K = parameters[i][2];
//        int TILE = parameters[i][3];
//
//        // declare matrix
//        int* A = (int*)malloc(sizeof(int) * M * N);
//        int* B = (int*)malloc(sizeof(int) * N * K);
//        int* C_BruteForce = (int*)malloc(sizeof(int) * M * K);  // using brute force matrix multiplication
//        int* C_Tile = (int*)malloc(sizeof(int) * M * K);     // using tiling to make cache hot
//        int* C_Tile_Thread = (int*)malloc(sizeof(int) * M * K);   // using tiling and multiple threads
//        memset(C_Tile, 0, M*K*sizeof(int));
//        memset(C_Tile_Thread, 0, M*K*sizeof(int));
//
//        generateMatrixUsingRandomNumbers(A, M, N, range);
//        generateMatrixUsingRandomNumbers(B, N, K, range);
//
//        // matrix AxB multiplication.
//        // Three Algorithms: algorithm 0: brute force; algorithm 1: using tiling; algorithm 2: tiling + multi_thread
//        long long exec_time_brute_force = getExecTime(A, B, C_BruteForce, M, N, K, TILE, 0);
//        long long exec_time_tile_only = getExecTime(A, B, C_Tile, M, N, K, TILE, 1);
//        long long exec_time_tile_threads = getExecTime(A, B, C_Tile_Thread, M, N, K, TILE, 2);
//
//        // will print "not equal" if any elelment is not equal for the three arrays
//        compareArrays(C_BruteForce, C_Tile, C_Tile_Thread, M, K);
//
//        printf("Dimensions of matrixes M=%d, N=%d, K=%d, tile_size=%d, threads=%d\n", M, N, K, TILE, THREAD_NUM);
//        printf("Matrix multiply exec time for brute force: %lld milliseconds\n", exec_time_brute_force);
//        printf("Matrix multiply exec time for tile method: %lld milliseconds\n", exec_time_tile_only);
//        printf("Matrix multiply exec time for tile+thread: %lld milliseconds\n\n", exec_time_tile_threads);
//    }
}

//// utility to print a matrix
//void print(int* matrix, int row, int col) {
//    for (int i = 0; i < row; i++) {
//        for (int j = 0; j < col; j++) {
//            printf("%*d ", WIDTH, matrix[i * col + j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//}

//void generateMatrixUsingRandomNumbers(int* matrix, int dim1, int dim2, int range) {
//    // Generate Matrix using random number in range [0,range-1] inclusive
//    for (int i = 0; i < dim1; ++i) {
//        for (int j = 0; j < dim2; ++j) {
//            matrix[i * dim2 + j] = rand() % range;
//        }
//    }
////    print(A, M, N);
//}

long long getExecTime(int* A, int* B, int* C, int M, int N, int K, int TILE, int algorithm) {
    auto start = std::chrono::high_resolution_clock::now();
    if (algorithm == 0) {
        bruteForceMultiply(A, B, C, M, N, K);
    }
    else if (algorithm == 1) {
        tiledMultiply(A, B, C, M, N, K, TILE);
    }
    else if (algorithm == 2) {
        tiledMultiply_multiThread(A, B, C, M, N, K, TILE);
    }
//    print(C, M, K);   // print results of matrix multiplication
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() / 1000; // in ms
//    printf("exec time for algorithm %d: %lld milliseconds\n", algorithm, duration);
    return duration;
}

void bruteForceMultiply(int *A, int *B, int *C, int M, int N, int K) {
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

/*
 * Since this method will udpate all over the result matrix, it will have race condition with the other worker threads
 * that calculate the "tile-multiple-part" of the matrix.  Since this method calculates the residual blocks which is
 * relatively small when matrices are large, it is safe and efficient to finish this work before starting other threads.
 */
void residualBlockThread(int* A, int* B, int* C, int M, int N, int K, int TILE) {
    // there may be residual sides along the bottom and right part of matrix A and/or B
    // first calculate indices for residual blocks
    // Note: this can't be declared at the beginning of this function and initialize them to 0. The reason is for the
    // first three loops above, the variable has to be initialize to zero for the 2nd and 3rd loop.
    int i = M - M % TILE;
    int j = N - N % TILE;
    int h = K - K % TILE;
    /*
     * First, for additional cols on the right side of A, there will be equal amount of rows at the bottom of B,
     * this is equal to an independent, separate matrix multiplication. We can calculate this part first.
     */
    for (int r1 = 0; r1 < M; r1++) {
        for (int inner1 = j; inner1 < N; inner1++) {
            for (int c1 = 0; c1 < K; c1++) {
                C[r1 * K + c1] += A[r1 * N + inner1] * B[inner1 * K + c1];
            }
        }
    }

    /*
     * second, for additional rows at the bottom of matrix A(exclude the part that overlaps with the residual
     * rectangle at the right
     *
     */
    for (int r2 = i; r2 < M; r2++) {
        for (int inner2 = 0; inner2 < j; inner2++) {
            for (int c2 = 0; c2 < K; c2++) {
                C[r2 * K + c2] += A[r2 * N + inner2] * B[inner2 * K + c2];
            }
        }
    }

    /*
     * third, for additional cols at the right of matrix B(exclude the part that overlaps with the residual
     * rectangle at the bottom
     *
     */
    for (int r3 = 0; r3 < i; r3++) {
        for (int inner3 = 0; inner3 < j; inner3++) {
            for (int c3 = h; c3 < K; c3++) {
                C[r3 * K + c3] += A[r3 * N + inner3] * B[inner3 * K + c3];
            }
        }
    }
}

/*
 * each thread will deal with some consecutive rows of the result matrix, thus work is divided.  Since the update of
 * the result matrix is row based, different thread access different rows, there is no data race condition. No mutex
 * needed.
 */
void tilingPartThreading(int* A, int* B, int* C, int M, int N, int K, int TILE, int thread_id) {
    int total_tiles_dir1 = K / TILE;
    int num_tile_per_thread = total_tiles_dir1 / (THREAD_NUM - 1);
    int start_idx = (thread_id - 1) * num_tile_per_thread * TILE;
    int end_idx = thread_id * num_tile_per_thread * TILE;
    if (thread_id == THREAD_NUM - 1) end_idx = M + 1 - TILE;

    for (int i = start_idx; i < end_idx; i += TILE) {
        for (int j = 0; j + TILE - 1 < N; j += TILE) {
            for (int h = 0; h + TILE - 1 < K; h += TILE) {  // be sure to start from zero for every upper loop
                // calculating a square tile unit
                for (int r = i; r < i + TILE; r++) {
                    for (int inner = j; inner < j + TILE; inner++) {
                        for (int c = h; c < h + TILE; c++) {
                            C[r * K + c] += A[r * N + inner] * B[inner * K + c];
                        }
                    }
                }
            }
        }
    }
}

void tiledMultiply_multiThread(int* A, int* B, int* C, int M, int N, int K, int TILE) {
    std::thread workers[THREAD_NUM];
    workers[0] = std::thread(residualBlockThread, std::ref(A), std::ref(B), std::ref(C), M, N, K, TILE);
    workers[0].join();  // this is required to avoid data race for larger matrices
    for (int i = 1; i < THREAD_NUM; i++) {
        workers[i] = std::thread(tilingPartThreading, std::ref(A), std::ref(B), std::ref(C), M, N, K, TILE, i);
    }

    for (int i = 1; i < THREAD_NUM; ++i) {
        workers[i].join();
    }
}

void tiledMultiply(int* A, int* B, int* C, int M, int N, int K, int TILE) {
    for (int i = 0; i + TILE - 1 < M; i += TILE) {
        for (int j = 0; j + TILE - 1 < N; j += TILE) {
            for (int h = 0; h + TILE - 1 < K; h += TILE) {  // be sure to start from zero for every upper loop
                // calculating a square tile unit
                for (int r = i; r < i + TILE; r++) {
                    for (int inner = j; inner < j + TILE; inner++) {
                        for (int c = h; c < h + TILE; c++) {
                            C[r * K + c] += A[r * N + inner] * B[inner * K + c];
                        }
                    }
                }
            }
        }
    }


    // there may be residual sides along the bottom and right part of matrix A and/or B
    // first calculate indices for residual blocks
    // Note: this can't be declared at the beginning of this function and initialize them to 0. The reason is for the
    // first three loops above, the variable has to be initialize to zero for the 2nd and 3rd loop.
    int i = M - M % TILE;
    int j = N - N % TILE;
    int h = K - K % TILE;
    /*
     * First, for additional cols on the right side of A, there will be equal amount of rows at the bottom of B,
     * this is equal to an independent, separate matrix multiplication. We can calculate this part first.
     */
    for (int r1 = 0; r1 < M; r1++) {
        for (int inner1 = j; inner1 < N; inner1++) {
            for (int c1 = 0; c1 < K; c1++) {
                C[r1 * K + c1] += A[r1 * N + inner1] * B[inner1 * K + c1];
            }
        }
    }

    /*
     * second, for additional rows at the bottom of matrix A(exclude the part that overlaps with the residual
     * rectangle at the right
     *
     */
    for (int r2 = i; r2 < M; r2++) {
        for (int inner2 = 0; inner2 < j; inner2++) {
            for (int c2 = 0; c2 < K; c2++) {
                C[r2 * K + c2] += A[r2 * N + inner2] * B[inner2 * K + c2];
            }
        }
    }

    /*
     * third, for additional cols at the right of matrix B(exclude the part that overlaps with the residual
     * rectangle at the bottom
     *
     */
    for (int r3 = 0; r3 < i; r3++) {
        for (int inner3 = 0; inner3 < j; inner3++) {
            for (int c3 = h; c3 < K; c3++) {
                C[r3 * K + c3] += A[r3 * N + inner3] * B[inner3 * K + c3];
            }
        }
    }
}

//void compareArrays(int* m1, int *m2, int* m3, int row, int col) {
//    for (int i = 0; i < row; i++) {
//        for (int j = 0; j < col; j++) {
//            if (m1[i * col + j] != m2[i * col + j] || m1[i * col + j] != m3[i * col + j]) {
//                printf(" ====== The three matrices calculated are NOT equal ======\n");
////                printf("i=%d, j=%d\n", i, j);
////                printf("m1=%d, m2=%d\n", m1[i * col + j], m2[i * col + j]);
//                return;
//            }
//        }
//    }
//    printf(" ====== The three matrices calculated are equal ======\n");
//}