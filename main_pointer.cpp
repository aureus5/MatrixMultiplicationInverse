//
// Created by GuoweiWu on 8/21/18.
//
#include <iostream>
#include <thread>
#include <chrono>

#define M 350
#define N 750
#define K 240
#define TILE 24
#define THREAD_NUM 8
#define WIDTH 3

void print(int* matrix, int row, int col);
void bruteForceMultiply(int* A, int* B, int* C);
void tiledMultiply(int* A, int* B, int* C, int tileLen);
void tiledMultiply_multiThread(int* A, int* B, int* C);
void residualBlockThread(int* A, int* B, int* C);
void tilingPartThreading(int* A, int* B, int* C, int thread_id);
void compareArrays(int* m1, int *m2, int row, int col);
std::mutex m;

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
//    tiledMultiply(A, B, C_Tiled, TILE);
    tiledMultiply_multiThread(A, B, C_Tiled);
    auto end_tile = std::chrono::high_resolution_clock::now();
    auto duration_tile = std::chrono::duration_cast<std::chrono::microseconds>(end_tile - start_tile);
//    print(C_Tiled, M, K);

    compareArrays(C_BruteForce, C_Tiled, M, K);

    printf("exec time for brute force: %lld milliseconds\n", duration_bf.count()/1000);
    printf("exec time for tile method: %lld milliseconds\n", duration_tile.count()/1000);
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
                            m.lock();
                            C[r * K + c] += A[r * N + inner] * B[inner * K + c];
                            m.unlock();
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

    for (int i = 1; i < THREAD_NUM; i++) {
        workers[i] = std::thread(tilingPartThreading, std::ref(A), std::ref(B), std::ref(C), i);
    }

    for (int i = 0; i < THREAD_NUM; ++i) {
        workers[i].join();
    }
}

void tiledMultiply(int* A, int* B, int* C, int tileLen) {
    for (int i = 0; i + tileLen - 1 < M; i += tileLen) {
        for (int j = 0; j + tileLen - 1 < N; j += tileLen) {
            for (int h = 0; h + tileLen - 1 < K; h += tileLen) {  // be sure to start from zero for every upper loop
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


    // there may be residual sides along the bottom and right part of matrix A and/or B
    // first calculate indices for residual blocks
    // Note: this can't be declared at the beginning of this function and initialize them to 0. The reason is for the
    // first three loops above, the variable has to be initialize to zero for the 2nd and 3rd loop.
    int i = M - M % tileLen;
    int j = N - N % tileLen;
    int h = K - K % tileLen;
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

void compareArrays(int* m1, int *m2, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (m1[i * col + j] != m2[i * col + j]) {
                printf(" ====== Not equal======\n");
                printf("i=%d, j=%d\n", i, j);
                printf("m1=%d, m2=%d\n", m1[i * col + j], m2[i * col + j]);
                return;
            }
        }
    }
    printf(" ====== equal======\n");
}