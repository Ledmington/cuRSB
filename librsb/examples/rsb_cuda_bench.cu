#include <stdio.h> /* printf() */
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "../rsb_cuda.h"

#define RSB_CHECK_ERROR(err)         \
    do                               \
    {                                \
        if (err != RSB_ERR_NO_ERROR) \
        {                            \
            rsb_perror(NULL, err);   \
            exit(EXIT_FAILURE);      \
        }                            \
    } while (0)

int main(const int argc, char *const argv[])
{
    srand(time(NULL));
    const rsb_blk_idx_t bs = RSB_DEFAULT_BLOCKING;
    const rsb_blk_idx_t brA = bs, bcA = bs;
    const RSB_DEFAULT_TYPE one = 1;
    const rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;

    const rsb_nnz_idx_t nnzA = 5; /* matrix nonzeroes count */
    const rsb_coo_idx_t nrA = 3;  /* matrix rows count */
    const rsb_coo_idx_t ncA = 3;  /* matrix columns count */

    /* nonzero row indices coordinates: */
    const rsb_coo_idx_t IA[] = {0, 0, 1, 1, 2};

    /* nonzero column indices coordinates: */
    const rsb_coo_idx_t JA[] = {0, 1, 1, 2, 0};

    const RSB_DEFAULT_TYPE VA[] = {9, 2, 11, 36, 13}; /* values of nonzeroes */

    RSB_DEFAULT_TYPE X[] = {0, 0, 0};          /* X vector's array */
    const RSB_DEFAULT_TYPE B[] = {42, 70, 16}; /* B vector's array */
    const size_t result_size = sizeof(X) / sizeof(X[0]);
    char ib[200];
    struct rsb_mtx_t *mtxAp = NULL; /* matrix structure pointer */

    cudaEvent_t start, end;
    float milliseconds = 0;
    cudaCheckError(cudaEventCreate(&start));
    cudaCheckError(cudaEventCreate(&end));

    rsb_err_t errval = RSB_ERR_NO_ERROR;

    printf("Hello, cuRSB!\n");
    printf("Initializing the library...\n");
    RSB_CHECK_ERROR(rsb_lib_init(RSB_NULL_INIT_OPTIONS));
    printf("Correctly initialized the library.\n");

    printf("Allocating matrix.\n");
    mtxAp = rsb_cuda_mtx_alloc_from_coo_const(
        VA, IA, JA, nnzA, typecode, nrA, ncA, brA, bcA,
        RSB_FLAG_NOFLAGS              /* default format will be chosen */
            | RSB_FLAG_DUPLICATES_SUM /* duplicates will be summed */
        ,
        &errval);
    if ((!mtxAp) || (errval != RSB_ERR_NO_ERROR))
    {
        printf("Error while allocating the matrix!\n");
        goto err;
    }
    printf("Correctly allocated a matrix.\n");

    printf("Summary information of the matrix:\n");
    /* print out the matrix summary information  */
    rsb_mtx_get_info_str(mtxAp, "RSB_MIF_MATRIX_INFO__TO__CHAR_P", ib, sizeof(ib));
    printf("%s", ib);
    printf("\n");

    printf("Computing SpMV.\n");
    cudaCheckError(cudaEventRecord(start));
    RSB_CHECK_ERROR(rsb_cuda_spmv(RSB_TRANSPOSITION_N, &one, mtxAp, B, 1, &one, X, 1));
    cudaCheckError(cudaDeviceSynchronize());
    cudaCheckError(cudaEventRecord(end));
    printf("Correctly performed a SpMV.\n");

    cudaEventElapsedTime(&milliseconds, start, end);
    printf("Elapsed time: %f seconds (%e)\n", milliseconds / 1000, milliseconds / 1000);

    printf("Printing result.\n");
    printf("[%f", X[0]);
    for (size_t i = 1; i < result_size; i++)
    {
        printf(", %f", X[i]);
    }
    printf("]\n");

    printf("Deallocating matrix.\n");
    rsb_cuda_mtx_free(mtxAp);
    printf("Correctly deallocated the matrix.\n");

    RSB_CHECK_ERROR(rsb_lib_exit(RSB_NULL_EXIT_OPTIONS));

    cudaCheckError(cudaEventDestroy(start));
    cudaCheckError(cudaEventDestroy(end));

    printf("Correctly finalized the library.\n");
    printf("Program terminating with no error.\n");
    return EXIT_SUCCESS;

err:
    rsb_perror(NULL, errval);
    printf("Program terminating with error.\n");
    return EXIT_FAILURE;
}
