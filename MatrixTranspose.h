//
// Created by GuoweiWu on 8/23/18.
//

#ifndef BRAINCORP_MATRIXTRANSPOSE_H
#define BRAINCORP_MATRIXTRANSPOSE_H

/*
 * brute force matrix transpose
 */
long long tileTransposeBruteForce(int *src, int *dest, int M, int N) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int row = 0; row < M; row++) {
        for (int col = 0; col < N; col++) {
            dest[col * M + row] = src[row * N + col];
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() / 1000; // in ms
    return duration;
}

/*
 * method to transpose matrix using tiling only. This generally results in <2 fold increase in efficiency
 */
long long tileTranspose(int *src, int *dest, int M, int N, int TILE) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int col = 0; col < N; col += TILE) {
        for (int row = 0; row < M; row += TILE) {
            for (int i = row; (i < row + TILE) && (i < M); i++) {
                for (int j = col; (j < col + TILE) && (j < N); j++) {
                    dest[j * M + i] = src[i * N + j];
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() / 1000; // in ms
    return duration;
}

/*
 * method for thread work. Work is divided by rows for dest matrix, thus there is no data race and no mutex needed
 */
void tilingThreadingUtil(int *src, int *dest, int M, int N, int TILE, int threads_num, int thread_id) {
    int total_tiles_dir1 = N / TILE;
    int num_tile_per_thread = total_tiles_dir1 / threads_num;
    int start_idx = thread_id * num_tile_per_thread * TILE;
    int end_idx = (thread_id + 1) * num_tile_per_thread * TILE;
    if (thread_id == threads_num - 1) end_idx = N;

    for (int col = start_idx; col < end_idx; col += TILE) {
        for (int row = 0; row < M; row += TILE) {
            for (int i = row; (i < row + TILE) && (i < M); i++) {
                for (int j = col; (j < col + TILE) && (j < N); j++) {
                    dest[j * M + i] = src[i * N + j];
                }

            }
        }
    }
};

/*
 * method to organize threads for dividing the work. Tiling+thread results in 3-10 fold increase in efficiency.
 */
long long tilePlusThreadTranspose(int *src, int *dest, int M, int N, int threads_num, int TILE) {
    auto start = std::chrono::high_resolution_clock::now();
    std::thread workers[threads_num];
    for (int i = 0; i < threads_num; i++) {
        workers[i] = std::thread(tilingThreadingUtil, std::ref(src), std::ref(dest), M, N, TILE, threads_num, i);
    }

    for (int i = 0; i < threads_num; ++i) {
        workers[i].join();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() / 1000; // in ms
    return duration;
}

/*
 * driver method for comparing different transpose methods.
 * Three methods are compared: 1. brute force, 2. using tiling only, 3. using tiles and multi-threads
 */
void transposeMatrix() {
    // initialize random generator
    srand(2018);
    int range = 100; // matrix elements are randoms numbers in the the range[0, range-1], inclusive

    int parameters[][4] = {{635,   783,   20},
                           {834,   895,   20},
                           {1376,  1487,  20},
                           {1687,  1589,  20},
                           {4376,  6487,  30},
                           {9687,  9589,  30},
                           {16387, 15489, 30}};

    printf("\n\nStarting of matrix transpose....\n\n");
    for (int i = 0; i < sizeof(parameters) / sizeof(parameters[0]); i++) {
        int M = parameters[i][0];
        int N = parameters[i][1];
        int TILE = parameters[i][2];

        // allocate memory for src and dest arrays
        int *src = (int *) malloc(sizeof(int) * M * N);
        int *dest_brute_force = (int *) malloc(sizeof(int) * N * M);
        int *dest_tile = (int *) malloc(sizeof(int) * N * M);
        int *dest_tile_thread = (int *) malloc(sizeof(int) * N * M);

        generateMatrixUsingRandomNumbers(src, M, N, range);

        long long exec_time_brute_force = tileTransposeBruteForce(src, dest_brute_force, M, N);
        long long exec_time_tiling = tileTranspose(src, dest_tile, M, N, TILE);
        int threads_num = 8;
        long long exec_time_tiling_thread = tilePlusThreadTranspose(src, dest_tile_thread, M, N, threads_num, TILE);

        // will print "not equal" if any elelment is not equal for the three arrays
        compareArrays(dest_brute_force, dest_tile, dest_tile_thread, N, M);

        printf("Matirx transpose. M=%d, N=%d, TILE=%d, absolute exec times are:\n", M, N, TILE);
        printf("brute force: %lld ms\ntiling only: %lld ms\ntiling_thread:%lld ms\n", exec_time_brute_force,
               exec_time_tiling, exec_time_tiling_thread);
        printf("fold of increase in efficiency: tiling/bf: %f, tiling_thread/bf: %f\n\n",
               exec_time_tiling == 0 ? -1 : (double) exec_time_brute_force / exec_time_tiling,
               exec_time_tiling_thread == 0 ? -1 : (double) exec_time_brute_force / exec_time_tiling_thread);
    }

    printf("\n\nEnd of matrix transpose....\n\n\n");
    for (int i = 0; i < 50; i++) {
        printf("=");
    }
    printf("\n\n\n");
}

#endif //BRAINCORP_MATRIXTRANSPOSE_H
