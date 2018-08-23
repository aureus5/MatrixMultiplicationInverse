//
// Created by GuoweiWu on 8/23/18.
//

#ifndef BRAINCORP_COMMONUTILS_H
#define BRAINCORP_COMMONUTILS_H
#define WIDTH 3

/*
 * utility to generate a matrix of specified dimensions using random numbers. The random generator has been
 * initialized earlier. Random numbers are in the range of [0, range-1], inclusive.
 */
void generateMatrixUsingRandomNumbers(int* matrix, int dim1, int dim2, int range) {
    // Generate Matrix using random number in range [0,range-1] inclusive
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            matrix[i * dim2 + j] = rand() % range;
        }
    }
//    print(A, M, N);
}


/*
 * utility to print a matrix
 */
void print(int* matrix, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%*d ", WIDTH, matrix[i * col + j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*
 * utility to compare if three arrays are identical. If any elelments is not the same among the three, print
 * "NOT equal", otherwise print "equal" at the end.
 */
void compareArrays(int* m1, int *m2, int* m3, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (m1[i * col + j] != m2[i * col + j] || m1[i * col + j] != m3[i * col + j]) {
                printf(" ====== The three matrices calculated are NOT equal ======\n");
//                printf("i=%d, j=%d\n", i, j);
//                printf("m1=%d, m2=%d\n", m1[i * col + j], m2[i * col + j]);
                return;
            }
        }
    }
    printf(" ====== The three matrices calculated are equal ======\n");
}

#endif //BRAINCORP_COMMONUTILS_H
