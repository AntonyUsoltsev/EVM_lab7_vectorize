#include <stdio.h>
#include <malloc.h>
#include <math.h>

#define N 4
#define M 10

void mult(float *B, float *A, float *R) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float summ = 0;
            for (int k = 0; k < N; k++) {
                summ += B[i * N + k] * A[j * N + k];
            }
            R[i * N + j] = summ;
        }
    }
}

void summ(float *A, float *B, float *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i * N + j] = A[i * N + j] + B[i * N + j];
        }
    }
}

void print_matr(float *A) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", A[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main() {
    float *A = calloc(N * N, sizeof(float));
    float *A_1 = calloc(N * N, sizeof(float));
    float *I = calloc(N * N, sizeof(float));
    float *B = calloc(N * N, sizeof(float));
    float *R = calloc(N * N, sizeof(float));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j)
                A[i * N + j] = 4.0f;
            else
                A[i * N + j] = 5.0f;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j)
                I[i * N + j] = 1.0f;
            else
                I[i * N + j] = 0.0f;
        }
    }

    float a_1 = -1;
    float a_inf = -1;
    float tmp = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            tmp += fabsf(A[j * N + i]);
        }
        if (tmp > a_1)
            a_1 = tmp;
        tmp = 0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            tmp += fabsf(A[i * N + j]);
        }
        if (tmp > a_inf)
            a_inf = tmp;
        tmp = 0;
    }
    float koef = a_1 * a_inf;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i * N + j] = A[j * N + i] / koef;
        }
    }
    //print_matr(B);

    mult(B, A, R);
    //print_matr(B);
    // print_matr(R);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = I[j * N + i] - R[i * N + j];
        }
    }
    // print_matr(R);

    float *cur = calloc(N * N, sizeof(float));
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            cur[i * N + j] = R[i * N + j];
//        }
//    }
    float *cur_2 = calloc(N * N, sizeof(float));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cur_2[i * N + j] = R[i * N + j];
        }
    }
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                cur[i * N + j] = cur_2[i * N + j];
            }
        }
        summ(I, cur, I);
        mult(cur, R, cur_2);

    }
    mult(I, B, A_1);
    print_matr(A_1);
    return 0;
}
