#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <xmmintrin.h>
#include <time.h>
#define N 8
#define M 100000
void mult(float *a, float *b, float *R) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            R[i * N + j] = 0;
    __m128 *xx = (__m128 *) b;
    __m128 *m128_result = (__m128 *) R;
    __m128 mult, tmp;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            mult = _mm_set1_ps(a[i * N + j]);
            for (int k = 0; k < N / 4; ++k) {
                tmp = _mm_mul_ps(mult, xx[N * j / 4 + k]);
                m128_result[N * i / 4 + k] = _mm_add_ps(m128_result[N * i / 4 + k], tmp);
            }
        }
    }
}
//void mult(float *B, float *A, float *R) {
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            float summ = 0;
//            for (int k = 0; k < N; k++) {
//                summ += B[i * N + k] * A[k * N + j];
//            }
//            R[i * N + j] = summ;
//        }
//    }
//}

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
    float *cur = calloc(N * N, sizeof(float));
    float *cur_2 = calloc(N * N, sizeof(float));
    clock_t beg = clock();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j)
                A[i * N + j] = 4.0f;
            else
                A[i * N + j] = 3.0f;
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
     //print_matr(R);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = I[i * N + j] - R[i * N + j];
        }
    }
    // print_matr(R);


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
   // print_matr(A_1);
    mult(A,A_1,cur);
    //print_matr(cur);
    clock_t end = clock();
    double time_used = ((double) (end - beg)) / CLOCKS_PER_SEC;
    printf("%lf\n", time_used);
    return 0;
}
