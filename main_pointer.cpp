//
// Created by GuoweiWu on 8/21/18.
//
#include <iostream>
#include <thread>
#include <chrono>

#define M 1350
#define N 1750
#define K 1240
#define TILE 24
#define THREAD_NUM 8
#define WIDTH 3

void print(int* matrix, int row, int col);
void generateMatrixUsingRandomNumbers(int* matrix, int dim1, int dim2, int range);
long long getExecTime(int* A, int* B, int* C, int algorithm);
void bruteForceMultiply(int* A, int* B, int* C);
void tiledMultiply(int* A, int* B, int* C);
void tiledMultiply_multiThread(int* A, int* B, int* C);
void residualBlockThread(int* A, int* B, int* C);
void tilingPartThreading(int* A, int* B, int* C, int thread_id);
void compareArrays(int* m1, int* m2, int* m3, int row, int col);

int main() {

    // initialize random generator
    srand(2018);

    int range = 10;

    // declare matrix
    int* A = (int*)malloc(sizeof(int) * M * N);
    int* B = (int*)malloc(sizeof(int) * N * K);
    int* C_BruteForce = (int*)malloc(sizeof(int) * M * K);  // using brute force matrix multiplication
    int* C_Tile = (int*)malloc(sizeof(int) * M * K);     // using tiling to make cache hot
    int* C_Tile_Thread = (int*)malloc(sizeof(int) * M * K);   // using tiling and multiple threads
    memset(C_Tile, 0, M*K*sizeof(int));
    memset(C_Tile_Thread, 0, M*K*sizeof(int));

    generateMatrixUsingRandomNumbers(A, M, N, range);
    generateMatrixUsingRandomNumbers(B, N, K, range);

    long long exec_time_brute_force = getExecTime(A, B, C_BruteForce, 0);
    long long exec_time_tile_only = getExecTime(A, B, C_Tile, 1);
    long long exec_time_tile_threads = getExecTime(A, B, C_Tile_Thread, 2);

    // will print "not equal" if any elelment is not equal for the three arrays
    compareArrays(C_BruteForce, C_Tile, C_Tile_Thread, M, K);

    printf("exec time for brute force: %lld milliseconds\n", exec_time_brute_force);
    printf("exec time for tile method: %lld milliseconds\n", exec_time_tile_only);
    printf("exec time for tile_thread method: %lld milliseconds\n", exec_time_tile_threads);
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

void generateMatrixUsingRandomNumbers(int* matrix, int dim1, int dim2, int range) {
    // Generate Matrix using random number in range [0,range-1] inclusive
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            matrix[i * dim2 + j] = rand() % range;
        }
    }
//    print(A, M, N);
}

long long getExecTime(int* A, int* B, int* C, int algorithm) {
    auto start = std::chrono::high_resolution_clock::now();
    if (algorithm == 0) {
        bruteForceMultiply(A, B, C);
    }
    else if (algorithm == 1) {
        tiledMultiply(A, B, C);
    }
    else if (algorithm == 2) {
        tiledMultiply_multiThread(A, B, C);
    }
//    print(C, M, K);   // print results of matrix multiplication
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_bf = (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() / 1000; // in ms
    return duration_bf;
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

void residualBlockThread(int* A, int* B, int* C) {
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

void tilingPartThreading(int* A, int* B, int* C, int thread_id) {
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

void tiledMultiply_multiThread(int* A, int* B, int* C) {
    std::thread workers[THREAD_NUM];
    workers[0] = std::thread(residualBlockThread, std::ref(A), std::ref(B), std::ref(C));
    workers[0].join();
    for (int i = 1; i < THREAD_NUM; i++) {
        workers[i] = std::thread(tilingPartThreading, std::ref(A), std::ref(B), std::ref(C), i);
    }

    for (int i = 1; i < THREAD_NUM; ++i) {
        workers[i].join();
    }
}

void tiledMultiply(int* A, int* B, int* C) {
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

void compareArrays(int* m1, int *m2, int* m3, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (m1[i * col + j] != m2[i * col + j] || m1[i * col + j] != m3[i * col + j]) {
                printf(" ====== Not equal ======\n");
//                printf("i=%d, j=%d\n", i, j);
//                printf("m1=%d, m2=%d\n", m1[i * col + j], m2[i * col + j]);
                return;
            }
        }
    }
    printf(" ====== equal ======\n");
}