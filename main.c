#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <xmmintrin.h>
#include <time.h>
#define N 8
#define M 10
void mult1(float *a, float *b, float *R) {
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
void mult2(const float *A, const float *B, float *R) {
    float *copy = calloc(N * N, sizeof(float));;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            copy[i * N + j] = B[j * N + i];
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float summ = 0;
            for (int k = 0; k < N; k++) {
                summ += A[i * N + k] * copy[j * N + i];
            }
            R[i * N + j] = summ;
        }
    }
    free(copy);
}

void mult(float *A, float *B, float *R) {
    float *copy = calloc(N * N, sizeof(float));;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            copy[i * N + j] = B[j * N + i];
        }
    }
    __m128 *m128_left = (__m128 *) A;
    __m128 *m128_right = (__m128 *) copy;
    __m128 *m128_result = (__m128 *) R;
    __m128 mult;

    __m128 tmp;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j]=0;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // mult = _mm_set1_ps(A[i * N + j]);
            __m128 s = _mm_setzero_ps();

            for (int k = 0; k < N / 4; k++) {
                tmp = _mm_mul_ps(m128_left[i * N + k], m128_right[j * N + k]);
                s = _mm_add_ps(s, tmp);
            }
            tmp = _mm_movehl_ps(tmp, s);
            s = _mm_add_ps(s, tmp);
            tmp = _mm_shuffle_ps(s, s, 1);
            s = _mm_add_ss(s, tmp);
            float sum;
            _mm_store_ss(&sum, s);
            R[i * N + j] = sum;
            // return sum;

        }
    }
}


void summ1(float *A, float *B, float *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i * N + j] = A[i * N + j] + B[i * N + j];
        }
    }
}
void summ(const float *A, const float *B, float *C) {
    __m128 *m128_result = (__m128 *) C;
    __m128 *m128_left = (__m128 *) A;
    __m128 *m128_right = (__m128 *) B;
    for (int i = 0; i < N * N / 4; ++i) {
        m128_result[i] = _mm_add_ps(m128_left[i], m128_right[i]);
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
    float trace;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i==j)
                trace+=A_1[i*N+j];
        }
    }
    double time_used = ((double) (end - beg)) / CLOCKS_PER_SEC;
    printf("Time: %lf\n", time_used);
    printf("Trace: %f\n", trace);
    print_matr(A_1);
    return 0;
}
////Vect Time: 34.513000
////Trace :0.000759

//// Classic time : 120
//// Trace :
