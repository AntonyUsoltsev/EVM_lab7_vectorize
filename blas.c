#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "mkl.h"

#define N 2048
#define M 10

void create_matrix(float *Matr, float value1, float value2) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (i == j)
                Matr[i * N + j] = value1;
            else
                Matr[i * N + j] = value2;
}

float find_max_str(const float *A) {
    float a_inf = -1;
    float tmp = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            tmp += fabsf(A[i * N + j]);
        if (tmp > a_inf)
            a_inf = tmp;
        tmp = 0;
    }
    return a_inf;
}

float find_max_column(const float *A) {
    float a_1 = -1;
    float tmp = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++)
            tmp += fabsf(A[j * N + i]);
        if (tmp > a_1)
            a_1 = tmp;
        tmp = 0;
    }
    return a_1;
}

void create_B_matrix(const float *A, float *B, float a_1, float a_inf) {
    float koef = a_1 * a_inf;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            B[i * N + j] = A[j * N + i] / koef;
}

void transposition(const float *A, float *A_t) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A_t[i * N + j] = A[j * N + i];
}

void mult (const float *A, const float *B, float *C) {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, 1.0, A, N, B, N, 0.0, C, N);
}

void summ(const float *A, const float *B, float *C) {
    for (int i = 0; i < N * N; i++)
        C[i] = A[i] + B[i];
}

void print_matr(const float *A) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%lf ", A[i * N + j]);
        printf("\n");
    }
    printf("\n");
}

float count_trace(const float *A) {
    float trace = 0 ;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (i == j)
                trace += A[i * N + j];
    return trace;
}

int main() {
    float *A = calloc(N * N, sizeof(float));
    float *A_1 = calloc(N * N, sizeof(float));
    float *A_1_t = calloc(N * N, sizeof(float));
    float *A_t = calloc(N * N, sizeof(float));
    float *R_t = calloc(N * N, sizeof(float));
    float *B_t = calloc(N * N, sizeof(float));
    float *I = calloc(N * N, sizeof(float));
    float *B = calloc(N * N, sizeof(float));
    float *R = calloc(N * N, sizeof(float));
    float *cur = calloc(N * N, sizeof(float));
    float *cur_2 = calloc(N * N, sizeof(float));

    clock_t beg = clock();

    create_matrix(A, 4.0f, 3.0f);
    create_matrix(I, 1.0f, 0.0f);

    float a_1 = find_max_column(A);
    float a_inf = find_max_str(A);

    create_B_matrix(A, B, a_1, a_inf);

    transposition(A, A_t);
    mult(B, A_t, R);

    for (int i = 0; i < N * N; i++)
        R[i] = I[i] - R[i];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            cur_2[i * N + j] = R[i * N + j];

    transposition(R, R_t);

    // cur == R^(m-1), cur_2 = R^m
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                cur[i * N + j] = cur_2[i * N + j];
        summ(I, cur, I);
        mult(cur, R_t, cur_2);
    }

    transposition(B, B_t);
    mult(I, B_t, A_1);

    clock_t end = clock();

    transposition(A_1, A_1_t);
    mult(A, A_1_t, cur);

    float trace = count_trace(cur);

    double time_used = ((double) (end - beg)) / CLOCKS_PER_SEC;
    printf("Time: %lf\n", time_used);
    printf("Trace :%f\n", trace);

    // print_matr(A_1);
    //print_matr(cur);

    free(A); free(A_1); free(A_t);
    free(A_1_t); free(B); free(B_t);
    free(I); free(R); free(R_t);
    free(cur); free(cur_2);
    return 0;
}

//free(A); free(A_1); free(A_t);
//free(A_1_t); free(B); free(B_t);
//free(I); free(R); free(R_t);
//free(cur); free(cur_2);