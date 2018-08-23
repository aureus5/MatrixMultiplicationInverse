//
// Created by GuoweiWu on 8/21/18.
//
#include <iostream>
#include <cstdlib>

#define M 8
#define N 7
#define K 10
#define WIDTH 3

// utility to print a matrix
void print_MN(int matrix[][N]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%*d ", WIDTH, matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// utility to print a matrix
void print_NK(int matrix[][K]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%*d ", WIDTH, matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}



int main() {
printf("%d\n\n",(int)1e5);
    int* a = (int*)malloc(3*sizeof(int));
    memset(a, 0, sizeof(a));
    printf("%d\n", a[0]);
    printf("%d\n", a[1]);
    printf("%d\n", a[2]);
    printf("%d\n", a[3]);
    printf("%d\n", a[5]);

    // initialize random generator
    srand(2018);

    // declare matrix
    int A[M][N];
    int B[N][K];

    // Generate Matrix using random number in range [0,99]
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = rand() % 100;
        }
    }

    print_MN(A);

    // Generate Matrix using random number in range [0,99]
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < K; ++j) {
            B[i][j] = rand() % 100;
        }
    }

    print_NK(B);


}