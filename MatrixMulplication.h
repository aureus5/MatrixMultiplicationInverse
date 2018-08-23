//
// Created by GuoweiWu on 8/23/18.
//

#ifndef BRAINCORP_MATRIXMULPLICATION_H
#define BRAINCORP_MATRIXMULPLICATION_H

#include <cstring>
#define THREAD_NUM 8

/*
 * brute force using three loops.
 */
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

/*
 * driver method for setting up worker threads to devide the work.
 * the first thread(thread[0]) is used to calculate the residual blocks. This has to finish first before other
 * worker threads start. The reason is that this first thread will be updating the same entries as some of the worker
 * threads, leading to data race.
 */
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

/*
 * multiply using tiling only.  Since most matrices are not multiples of tiles, there will be residual blocks.
 * Here I am dividing the work into four parts: the computation heavy one is the portion that are multiples of
 * tiles(for large matrices). The rest will update the result matrix with necessary additions.
 */
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

/*
 * Utility method to obtain exec time for the three multiply algorithms:
 * 0. brute force, 1. tiling only, 2. tiling+multithread.
 */
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

/*
 * driver method to compare three multiply algorithms: 0. brute force, 1. tiling only, 2. tiling+multithread.
 */
void matrixMulplication() {
    // initialize random generator
    srand(2018);
    int range = 100; // matrix elements are randoms numbers in the the range[0, range-1], inclusive

    int parameters[][4] = {{135,  175,  124,  20},
                           {374,  567,  421,  20},
                           {635,  783,  323,  20},
                           {834,  895,  798,  20},
                           {1034,  1195,  498,  20},
                           {1376, 1487, 1533, 30},
                           {1687, 1589, 1863, 30}};

    printf("\n\nStarting of matrix multiplication....\n\n");

    for (int i = 0; i < sizeof(parameters) / sizeof(parameters[0]); i++) {
        int M = parameters[i][0];
        int N = parameters[i][1];
        int K = parameters[i][2];
        int TILE = parameters[i][3];

        // declare matrix
        int* A = (int*)malloc(sizeof(int) * M * N);
        int* B = (int*)malloc(sizeof(int) * N * K);
        int* C_BruteForce = (int*)malloc(sizeof(int) * M * K);  // using brute force matrix multiplication
        int* C_Tile = (int*)malloc(sizeof(int) * M * K);     // using tiling to make cache hot
        int* C_Tile_Thread = (int*)malloc(sizeof(int) * M * K);   // using tiling and multiple threads
        std::memset(C_Tile, 0, M*K*sizeof(int));
        std::memset(C_Tile_Thread, 0, M*K*sizeof(int));

        generateMatrixUsingRandomNumbers(A, M, N, range);
        generateMatrixUsingRandomNumbers(B, N, K, range);

        // matrix AxB multiplication.
        // Three Algorithms: algorithm 0: brute force; algorithm 1: using tiling; algorithm 2: tiling + multi_thread
        long long exec_time_brute_force = getExecTime(A, B, C_BruteForce, M, N, K, TILE, 0);
        long long exec_time_tile_only = getExecTime(A, B, C_Tile, M, N, K, TILE, 1);
        long long exec_time_tile_threads = getExecTime(A, B, C_Tile_Thread, M, N, K, TILE, 2);

        // will print "not equal" if any elelment is not equal for the three arrays
        compareArrays(C_BruteForce, C_Tile, C_Tile_Thread, M, K);

        printf("Matrix multiplication. M=%d, N=%d, K=%d, tile_size=%d, threads=%d\n", M, N, K, TILE, THREAD_NUM);
        printf("Absolute matrix multiplication exec times are:\n");
        printf("brute force: %lld ms\ntiling only: %lld ms\ntiling_thread:%lld ms\n", exec_time_brute_force,
               exec_time_tile_only, exec_time_tile_threads);
        printf("fold of increase in efficiency: tiling/bf: %f, tiling_thread/bf: %f\n\n",
               exec_time_tile_only == 0 ? -1 : (double) exec_time_brute_force / exec_time_tile_only,
               exec_time_tile_threads == 0 ? -1 : (double) exec_time_brute_force / exec_time_tile_threads);

    }
    printf("\n\nEnd of matrix transpose....\n\n\n");
}

#endif //BRAINCORP_MATRIXMULPLICATION_H
