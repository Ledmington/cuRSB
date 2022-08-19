#ifndef RSB_CUDA_H_INCLUDED
#define RSB_CUDA_H_INCLUDED

#include "rsb.h"

#define cudaCheckError(cudaAPICall)                                                                                                              \
    {                                                                                                                                            \
        cudaError_t my_error_id = cudaAPICall;                                                                                                   \
        if (my_error_id != cudaSuccess)                                                                                                          \
        {                                                                                                                                        \
            /* __func__ is c99 and c++11 compatible but on older versions it may be necessary to use __FUNCTION__ instead */                     \
            printf("CUDA error at \"%s\":%s:%d = %d\n-> %s\n", __FILE__, __func__, __LINE__, (int)my_error_id, cudaGetErrorString(my_error_id)); \
            /*resetAllDevices();*/                                                                                                               \
            exit(EXIT_FAILURE);                                                                                                                  \
        }                                                                                                                                        \
    }

#define RSB_INTERFACE_RETURN_ERR(ERRVAL) RSB_INTERFACE_ENDCMD RSB_DO_ERR_RETURN_INTERFACE(ERRVAL)

// Memory management functions
struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flagsA, rsb_err_t *errvalp);
rsb_err_t rsb_cuda_mtx_alloc_from_coo_end(struct rsb_mtx_t **mtxApp);
struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_csr_const(const void *VA, const rsb_coo_idx_t *RP, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t *errvalp);
struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_csc_const(const void *VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *CP, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t *errvalp);
struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_csr_inplace(void *VA, rsb_nnz_idx_t *RP, rsb_coo_idx_t *JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t *errvalp);
struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t *errvalp);
struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_coo_inplace(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t *errvalp);
rsb_err_t rsb_cuda_mtx_clone(struct rsb_mtx_t **mtxBpp, rsb_type_t typecode, rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t *mtxAp, rsb_flags_t flags);
struct rsb_mtx_t *rsb_cuda_mtx_free(struct rsb_mtx_t *mtxAp);

// Computation functions
rsb_err_t rsb_cuda_spmv(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t *mtxAp, const void *Xp, rsb_coo_idx_t incX, const void *betap, void *Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_cuda_spmm(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t *mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void *Bp, rsb_nnz_idx_t ldB, const void *betap, void *Cp, rsb_nnz_idx_t ldC);
rsb_err_t rsb_cuda_spsv(rsb_trans_t transT, const void *alphap, const struct rsb_mtx_t *mtxTp, const void *Xp, rsb_coo_idx_t incX, void *Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_cuda_spsm(rsb_trans_t transT, const void *alphap, const struct rsb_mtx_t *mtxTp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void *betap, const void *Bp, rsb_nnz_idx_t ldB, void *Cp, rsb_nnz_idx_t ldC);

#endif // RSB_CUDA_H_INCLUDED