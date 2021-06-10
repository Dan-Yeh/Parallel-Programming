#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <immintrin.h>

#pragma GCC optimize ("-Ofast")

void dgemm(float alpha, const float *a, const float *b, float beta, float *c, int paddedSize) {
    int AVX_SIZE = MATRIX_SIZE - MATRIX_SIZE % 8;
    float tmp_sum = 0;

    #pragma omp parallel for firstprivate(AVX_SIZE, tmp_sum) shared(a, b, c) schedule(dynamic, 16) num_threads(32)
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            __m256 partial_sum = _mm256_set1_ps(0);
            for (int k = 0; k < AVX_SIZE; k+=8) {
                __m256 partial_a = _mm256_load_ps(a + i * paddedSize + k);
                __m256 partial_b = _mm256_load_ps(b + j * paddedSize + k);
                partial_a = _mm256_mul_ps(partial_a, partial_b);
                partial_sum = _mm256_add_ps(partial_sum, partial_a);
            }
            // _mm256_storeu_ps(partial_sum_array, partial_sum);
            partial_sum = _mm256_hadd_ps(partial_sum, partial_sum);
            // 0+1 0+1 2+3 2+3 4+5 4+5 6+7 6+7
            partial_sum = _mm256_hadd_ps(partial_sum, partial_sum);
            // 0+1+2+3 0+1+2+3 ... 4+5+6+7 ...
            __m256 permutePartialSum = _mm256_permute2f128_ps(partial_sum, partial_sum, 1);
            partial_sum = _mm256_add_ps(partial_sum, permutePartialSum);
            tmp_sum = _mm256_cvtss_f32(partial_sum);

            tmp_sum += alpha * a[i * MATRIX_SIZE + 2016] * b[j * MATRIX_SIZE + 2016];
            c[i * MATRIX_SIZE + j] = alpha * tmp_sum;
            tmp_sum = 0;
        }
    }
}


void generateProblemFromInput(float& alpha, float* a, float* b, float& beta, float* c, int paddedSize) {
    unsigned int seed = 0;
    std::cout << "READY" << std::endl;
    std::cin >> seed;

    std::cerr << "Using seed " << seed << std::endl;
    if (seed == 0) {
        std::cerr << "Warning: default value 0 used as seed." << std::endl;
    }

    std::mt19937 random(seed);
    std::uniform_real_distribution<float> distribution(-5, 5);

    /* initialisation */
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            *(a + i * paddedSize + j) = distribution(random);
            *(b + i * paddedSize + j) = distribution(random);
        }
    }
    memset(c, 0, MEM_SIZE);

    alpha = distribution(random);
    beta = distribution(random);
}


int main(int, char **) {
    float alpha, beta;

    // mem allocations
    int paddedSize = MATRIX_SIZE + ((32 - MATRIX_SIZE*sizeof(float)%32)%32)/sizeof(float);
                                    // (32 - (2017 * 4) %32) / 4 = 7
    int padded_mem_size = paddedSize * MATRIX_SIZE * sizeof(float);
    int mem_size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
    auto a = (float *) aligned_alloc(32, padded_mem_size);
    auto b = (float *) aligned_alloc(32, padded_mem_size);
    auto c = (float *) malloc(mem_size);

    // check if allocated
    if (nullptr == a || nullptr == b || nullptr == c) {
        printf("Memory allocation failed\n");
        if (nullptr != a) free(a);
        if (nullptr != b) free(b);
        if (nullptr != c) free(c);
        return 0;
    }

    generateProblemFromInput(alpha, a, b, beta, c, paddedSize);

    std::cerr << "Launching dgemm step." << std::endl;
    // matrix-multiplication
    dgemm(alpha, a, b, beta, c, paddedSize);

    outputSolution(c);

    free(a);
    free(b);
    free(c);
    return 0;
}
