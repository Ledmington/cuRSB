#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> // fabs
#include <assert.h>
#include <omp.h> // omp_get_wtime

#define EPS 1e-12
#define max(a,b) (a>b)?a:b

typedef struct {
    unsigned int rows;
    unsigned int columns;
    float *v;
} matrix;

typedef struct {
    unsigned int len;
    unsigned int *r_index;
    unsigned int *c_index;
    float *values;
} coo_matrix;

typedef struct {
    unsigned int len;
    unsigned int len_rows;
    unsigned int *r_index;
    unsigned int *c_index;
    float *values;
} csr_matrix;

float randab(const float a, const float b) {
    return (float)rand() / (float)RAND_MAX * (b-a) + a;
}

unsigned int get_matrix_size(const matrix *m) {
    return sizeof(m->rows) + sizeof(m->columns) + sizeof(m->v[0]) * m->rows * m->columns;
}

unsigned int get_coo_matrix_size(const coo_matrix *m) {
    return sizeof(m->len) +
            sizeof(m->c_index[0]) * m->len +
            sizeof(m->r_index[0]) * m->len +
            sizeof(m->values[0])  * m->len;
}

unsigned int get_csr_matrix_size(const csr_matrix *m) {
    return sizeof(m->len) +
            sizeof(m->len_rows) +
            sizeof(m->r_index[0]) * m->len_rows +
            sizeof(m->c_index[0]) * m->len +
            sizeof(m->values[0]) * m->len;
}

matrix* init_matrix(matrix *m, const unsigned int r, const unsigned int c) {
    m->rows = r;
    m->columns = c;
    m->v = (float*) calloc(r * c, sizeof(float));
    assert(m->v);
    const unsigned int n = max(r, c);
    for (unsigned int i=0; i<n; i++) {
        const unsigned int r_idx = rand() % r;
        const unsigned int c_idx = rand() % c;
        m->v[r_idx*m->columns+c_idx] = randab(0.0f, 10.0f);
    }
    return m;
}

coo_matrix* convert_to_coo(const matrix *m, coo_matrix *coo_m) {
    // Counting non-zero values
    coo_m->len = 0;
    for(unsigned int r=0; r<m->rows; r++) {
        for(unsigned int c=0; c<m->columns; c++) {
            if(fabs(m->v[r*m->columns+c]) > EPS) {
                coo_m->len++;
            }
        }
    }

    coo_m->r_index = (unsigned int*) malloc(coo_m->len * sizeof(unsigned int));
    assert(coo_m->r_index);
    coo_m->c_index = (unsigned int*) malloc(coo_m->len * sizeof(unsigned int));
    assert(coo_m->c_index);
    coo_m->values = (float*) malloc(coo_m->len * sizeof(float));
    assert(coo_m->values);

    // Filling the vectors
    unsigned int idx = 0;
    for(unsigned int r=0; r<m->rows; r++) {
        for(unsigned int c=0; c<m->columns; c++) {
            if(fabs(m->v[r*m->columns+c]) > EPS) {
                coo_m->r_index[idx] = r;
                coo_m->c_index[idx] = c;
                coo_m->values[idx] = m->v[r*m->columns+c];
                idx++;
            }
        }
    }

    return coo_m;
}

csr_matrix* convert_to_csr(const matrix *m, csr_matrix *csr_m) {
    // Counting non-zero values
    csr_m->len = 0;
    csr_m->len_rows = 0;
    int last_row_idx = -1;
    for(unsigned int r=0; r<m->rows; r++) {
        for(unsigned int c=0; c<m->columns; c++) {
            if(fabs(m->v[r*m->columns+c]) > EPS) {
                csr_m->len++;
                if (r != last_row_idx) {
                    csr_m->len_rows++;
                    last_row_idx = r;
                }
            }
        }
    }

    csr_m->r_index = (unsigned int*) malloc(csr_m->len_rows * sizeof(unsigned int));
    assert(csr_m->r_index);
    csr_m->c_index = (unsigned int*) malloc(csr_m->len * sizeof(unsigned int));
    assert(csr_m->c_index);
    csr_m->values = (float*) malloc(csr_m->len * sizeof(float));
    assert(csr_m->values);

    // Filling the vectors
    unsigned int idx = 0;
    unsigned int r_idx = 0;
    last_row_idx = -1;
    for(unsigned int r=0; r<m->rows; r++) {
        for(unsigned int c=0; c<m->columns; c++) {
            if(fabs(m->v[r*m->columns+c]) > EPS) {
                if(r != last_row_idx) {
                    csr_m->r_index[r_idx] = r;
                    r_idx++;
                    last_row_idx = r;
                }
                csr_m->c_index[idx] = c;
                csr_m->values[idx] = m->v[r*m->columns+c];
                idx++;
            }
        }
    }

    return csr_m;
}

void print_matrix(const matrix *m) {
    printf("%u x %u\n", m->rows, m->columns);
    for(unsigned int r=0; r<m->rows; r++) {
        printf("|");
        for(unsigned int c=0; c<m->columns; c++) {
            printf(" %.1f", m->v[r*m->columns+c]);
        }
        printf(" |\n");
    }
}

void print_coo_matrix(const coo_matrix *m) {
    printf("COO length: %u\n", m->len);
    printf("Row_idx: ");
    for(unsigned int i=0; i<m->len; i++) {
        printf("%3u ", m->r_index[i]);
    }
    printf("\nCol_idx: ");
    for(unsigned int i=0; i<m->len; i++) {
        printf("%3u ", m->c_index[i]);
    }
    printf("\nValues:  ");
    for(unsigned int i=0; i<m->len; i++) {
        printf("%.1f ", m->values[i]);
    }
    printf("\n");
}

void print_csr_matrix(const csr_matrix *m) {
    printf("CSR length: %u\n", m->len);
    printf("CSR rows length: %u\n", m->len_rows);
    printf("Row_idx: ");
    for(unsigned int i=0; i<m->len_rows; i++) {
        printf("%3u ", m->r_index[i]);
    }
    printf("\nCol_idx: ");
    for(unsigned int i=0; i<m->len; i++) {
        printf("%3u ", m->c_index[i]);
    }
    printf("\nValues:  ");
    for(unsigned int i=0; i<m->len; i++) {
        printf("%.1f ", m->values[i]);
    }
    printf("\n");
}

double benchmark_self_multiplication_normal(const matrix *m) {
    double t = omp_get_wtime();
    // TODO
    return omp_get_wtime() - t;
}

double benchmark_self_multiplication_coo(const coo_matrix *coo_m, const unsigned int rows, const unsigned int columns) {
    double t = omp_get_wtime();
    // TODO
    return omp_get_wtime() - t;
}

double benchmark_self_multiplication_csr(const csr_matrix *csr_m, const unsigned int rows, const unsigned int columns) {
    double t = omp_get_wtime();
    // TODO
    return omp_get_wtime() - t;
}

void destroy_matrix(const matrix *m) {
    free(m->v);
}

void destroy_coo_matrix(const coo_matrix *m) {
    free(m->c_index);
    free(m->r_index);
    free(m->values);
}

void destroy_csr_matrix(const csr_matrix *m) {
    free(m->c_index);
    free(m->r_index);
    free(m->values);
}

int main(const int argc, const char *argv[]) {
    srand(time(NULL));

    unsigned int rows = 12;
    unsigned int columns = 15;

    if(argc == 2 || argc > 3) {
        fprintf(stderr, "ERROR: wrong number of parameters\nUsage: %s [r c]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if(argc == 3) {
        rows = (unsigned int) strtoul(argv[1], (char**)NULL, 10);
        columns = (unsigned int) strtoul(argv[2], (char**)NULL, 10);
    }

    matrix m;
    coo_matrix coo_m;
    csr_matrix csr_m;

    init_matrix(&m, rows, columns);
    if(rows * columns < 100) print_matrix(&m);
    printf("Normal: %u (%.3e) bytes\n\n", get_matrix_size(&m), (float)get_matrix_size(&m));

    convert_to_coo(&m, &coo_m);
    if(rows * columns < 100) print_coo_matrix(&coo_m);
    printf("COO: %u (%.3e) bytes\n\n", get_coo_matrix_size(&coo_m), (float)get_coo_matrix_size(&coo_m));

    convert_to_csr(&m, &csr_m);
    if(rows * columns < 100) print_csr_matrix(&csr_m);
    printf("CSR: %u (%.3e) bytes\n\n", get_csr_matrix_size(&csr_m), (float)get_csr_matrix_size(&csr_m));

    printf("\n\n ### Benchmarks ###\n");
    printf("Self-multiplication\n");
    printf("Normal: %.3f seconds\n", benchmark_self_multiplication_normal(&m));
    printf("COO:    %.3f seconds\n", benchmark_random_access_coo(&coo_m, m.rows, m.columns));
    printf("CSR:    %.3f seconds\n", benchmark_random_access_csr(&csr_m, m.rows, m.columns));

    destroy_matrix(&m);
    destroy_coo_matrix(&coo_m);
    destroy_csr_matrix(&csr_m);

    return EXIT_SUCCESS;
}