//
// Created by GuoweiWu on 8/23/18.
//

#ifndef BRAINCORP_MATRIXTRANSPOSE_H
#define BRAINCORP_MATRIXTRANSPOSE_H

long long tileTransposeBruteForce(int* src, int* dest, int M, int N) {
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

long long tileTranspose(int* src, int* dest, int M, int N, int TILE) {
    auto start = std::chrono::high_resolution_clock::now();
    int diag = 0;
    int temp = 0;
    for (int col = 0; col < N; col += TILE) { /* same as above, iterating over columns */
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

void transposeMatrix(int M, int N, int TILE) {
    // initialize random generator
    srand(2018);
    int range = 100; // matrix elements are randoms numbers in the the range[0, range-1], inclusive

    int* src = (int*)malloc(sizeof(int) * M * N);
    int* dest_brute_force = (int*)malloc(sizeof(int) * N * M);
    int* dest_tile = (int*)malloc(sizeof(int) * N * M);
    int* dest_tile_thread = (int*)malloc(sizeof(int) * N * M);

    generateMatrixUsingRandomNumbers(src, M, N, range);
//    printf("before transpose: \n");
//    print(src, M, N);

    long long exec_time_brute_force = tileTransposeBruteForce(src, dest_brute_force, M, N);
//    printf("brute force: after transpose: \n");
//    print(dest_brute_force, N, M);

    long long exec_time_tiling = tileTranspose(src, dest_tile, M, N, TILE);
//    printf("tiling: after transpose: \n");
//    print(dest_tile, N, M);

    long long exec_time_tiling_thread = tileTranspose(src, dest_tile_thread, M, N, TILE);
//    printf("tiling: after transpose: \n");
//    print(dest_tile, N, M);


    printf("%lld\n%lld\n%lld\n", exec_time_brute_force, exec_time_tiling, exec_time_tiling_thread);

    // will print "not equal" if any elelment is not equal for the three arrays
    compareArrays(dest_brute_force, dest_tile, dest_tile_thread, N, M);
}

#endif //BRAINCORP_MATRIXTRANSPOSE_H
