#include <assert.h>

#define RSB_LOCK_H_INCLUDED // TODO: change or remove

#include "rsb_cuda.h"
#include "rsb_internals.h"
#include "rsb_do.h"
#include "rsb_common.h"
//#include "rsb_lock.h" // TODO remove if unused

#define min(A, B) ((A) < (B) ? (A) : (B))

__device__ inline rsb_bool_t rsb__is_equal(const void *x, rsb_type_t type, const void *y)
{
    switch (type)
    {
    case RSB_NUMERICAL_TYPE_DOUBLE:
        return (*(double *)x) == (*(double *)y);
        break;
    case RSB_NUMERICAL_TYPE_FLOAT:
        return (*(float *)x) == (*(float *)y);
        break;
    case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX:
        return (*(float _Complex *)x) == (*(float _Complex *)y);
        break;
    case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX:
        return (*(double _Complex *)x) == (*(double _Complex *)y);
        break;
    default:
        return 0;
        break;
    }
}

__device__ inline rsb_bool_t rsb__is_element_one(const void *x, rsb_type_t type)
{
    int y = 1;
    return rsb__is_equal(x, type, &y);
}

__device__ inline rsb_bool_t rsb__is_element_zero(const void *x, rsb_type_t type)
{
    int y = 0;
    return rsb__is_equal(x, type, &y);
}

__device__ inline rsb_bool_t rsb__is_element_minus_one(const void *x, rsb_type_t type)
{
    int y = -1;
    return rsb__is_equal(x, type, &y);
}

__device__ rsb_err_t rsb_cuda__do_spmv_non_recursive(const struct rsb_mtx_t *mtxAp, const void *x, void *y, const void *alphap, const void *betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
{
    rsb_err_t errval = RSB_ERR_NO_ERROR;

    const rsb_bool_t nostride = (incx == 1 && incy == 1) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;
    const rsb_bool_t should_scale_y = (betap && !rsb__is_element_one(betap, mtxAp->typecode)) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;
    const rsb_bool_t use_alpha_one = (!alphap || rsb__is_element_one(alphap, mtxAp->typecode)) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;
    const rsb_bool_t use_y_zeroing_kernel = (should_scale_y && rsb__is_element_zero(betap, mtxAp->typecode) && nostride && use_alpha_one) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;

    /*const size_t outtot=0,rhstot=0;
    const size_t outnri=0,rhsnri=0;
    const rsb_int_t nrhs=1;*/

    const size_t lenx = (mtxAp->el_size * rhsnri);
    const size_t leny = (mtxAp->el_size * outnri);
    rsb_int_t nrhsi = 0;

    for (nrhsi = 0; nrhsi < nrhs; ++nrhsi)
    {
        void *out = ((rsb_byte_t *)y) + (leny * nrhsi);
        const void *rhs = ((const rsb_byte_t *)x) + (lenx * nrhsi);

        if (should_scale_y && !use_y_zeroing_kernel)
        {
            RSB_CBLAS_X_SCAL_SPMV(mtxAp->typecode, rsb__do_get_rows_of(mtxAp, transA), betap, out, incy);
        }
        /* no beta specified counts as beta=1, and so no scaling is needed */

        if (use_alpha_one)
        {
            /* no alpha specified counts as alpha=1 */
            if (nostride)
            {
                if (use_y_zeroing_kernel)
                /* y <- a * x  */
                {
                    RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_uauz(mtxAp, rhs, out, transA));
                }
                else
                /* y <- y + a * x  */
                {
                    RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_uaua(mtxAp, rhs, out, transA));
                }
            }
            else
            /* y <- a * x  , with stride */
            {
                RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_sasa(mtxAp, rhs, out, incx, incy, transA));
            }
        }
        else
        {
            if (nostride)
            {
                /* y <- - a * x  */
                if (rsb__is_element_minus_one(alphap, mtxAp->typecode))
                {
                    RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_unua(mtxAp, rhs, out, transA));
                }
                /* y <- alpha * a * x  */
                else
                {
                    RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_uxua(mtxAp, rhs, out, alphap, transA));
                }
            }
            else
            /* y <- alpha * a * x  , with stride */
            {
                RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_sxsa(mtxAp, rhs, out, alphap, incx, incy, transA));
            }
        }
    }

    RSB_DO_ERR_RETURN(errval)
}

__global__ void rsb_cuda__spmv_kernel(const struct rsb_mtx_t *mtxAp, const void *x, void *y, const void *alphap, const void *betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPMV_ARGS)
{
    // const unsigned int thread_id_x = threadIdx.x;
    // const unsigned int thread_id_y = threadIdx.y;
    // const unsigned int block_id_x = blockIdx.x;
    // const unsigned int block_id_y = blockIdx.y;
    // const unsigned int block_dim_x = blockDim.x;
    // const unsigned int block_dim_y = blockDim.y;
    // const unsigned int grid_dim_x = gridDim.x;
    // const unsigned int grid_dim_y = gridDim.y;

    const unsigned int global_thread_id = blockIdx.x * blockDim.x + threadIdx.x;

    // printf("Thread (%u, %u) of (%u, %u) inside block (%u, %u) of (%u, %u)\n", thread_id_x, thread_id_y, block_dim_x, block_dim_y, block_id_x, block_id_y, grid_dim_x, grid_dim_y);
    printf("Thread %u (%d, %d, %d) alive\n", global_thread_id, threadIdx.x, blockDim.x, blockIdx.x);

    if (global_thread_id >= mtxAp->nnz)
    {
        printf("Thread %d dying\n", global_thread_id);
        return;
    }

    rsb_submatrix_idx_t n = 0;
    // rsb_submatrix_idx_t dm = 0; // TODO remove if unused/useless

    for (n = 0; n < mtxAp->all_leaf_matrices_n; ++n)
    {
        const struct rsb_mtx_t *const submatrix = mtxAp->all_leaf_matrices[n].mtxlp;
        rsb_bool_t gomv = RSB_BOOL_FALSE; // TODO remove if unused/useless

        char *const ov = (char *const)y;
        const rsb_coo_idx_t oincy = incy;

        // (CRITICAL SECTION)
        // {
        //     const rsb_coo_idx_t roff = submatrix->broff;
        //     const rsb_coo_idx_t nr = RSB_MTX_EFF_R(submatrix);
        //     const rsb_coo_idx_t coff = submatrix->bcoff;
        //     const rsb_coo_idx_t nc = RSB_MTX_EFF_C(submatrix);

        //     gomv = (rsb_do_spmv_lock_get(&lock, th_id, roff, nr, coff, nc, n, transA, &ov, &oincy) == RSB_BOOL_TRUE);
        //     if (gomv == RSB_BOOL_TRUE)
        //     {
        //         RSB_SPMV_VS_MARK_PRE(n);
        //     }
        // }

        // if (gomv == RSB_BOOL_TRUE)
        // {
        const size_t scoff = submatrix->coff - mtxAp->coff;
        const size_t sroff = submatrix->roff - mtxAp->roff;
        const char *const offx = ((const char *)x) + (el_size * scoff) * incx;
        char *const offy = ((char *)ov) + (el_size * sroff) * oincy;

        assert(scoff >= 0);
        assert(sroff >= 0);
        RSB_DO_ERROR_CUMULATE(errval, rsb_cuda__do_spmv_non_recursive(submatrix, offx, offy, alphap, NULL, incx, oincy, transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
        // }
    }

    // TODO return error somehow

    // #pragma omp parallel default(none)         \
//     shared(lock, all_leaf_matrices, mtxAp) \
//         reduction(|                        \
//                   : errval)
    //     {
    //         // const size_t want_threads = rsb_global_session_handle.rsb_want_threads;
    //         //  const size_t max_threads = min(all_leaf_matrices_n, want_threads);
    //         //  const rsb_thr_t th_id = omp_get_thread_num();
    //         rsb_submatrix_idx_t n = 0;
    //         rsb_submatrix_idx_t dm = 0;

    //         // if (RSB_UNLIKELY(th_id >= max_threads))
    //         // if (th_id >= max_threads)
    //         // {
    //         //     goto skip;
    //         // }
    //     again:
    //         // for (n = 0; RSB_LIKELY(n < all_leaf_matrices_n); ++n) {
    //         for (n = 0; n < all_leaf_matrices_n; ++n)
    //         {
    //             const struct rsb_mtx_t *const submatrix = all_leaf_matrices[n].mtxlp;
    //             rsb_bool_t gomv = RSB_BOOL_FALSE;

    //             char *const ov = (char *const)y;
    //             const rsb_coo_idx_t oincy = incy;

    // #pragma omp critical(rsb_spmv_crs)
    //             {
    //                 const rsb_coo_idx_t roff = submatrix->broff;
    //                 const rsb_coo_idx_t nr = RSB_MTX_EFF_R(submatrix);
    //                 const rsb_coo_idx_t coff = submatrix->bcoff;
    //                 const rsb_coo_idx_t nc = RSB_MTX_EFF_C(submatrix);

    //                 gomv = (rsb_do_spmv_lock_get(&lock, th_id, roff, nr, coff, nc, n, transA, &ov, &oincy) == RSB_BOOL_TRUE);
    //                 if (gomv == RSB_BOOL_TRUE)
    //                 {
    //                     RSB_SPMV_VS_MARK_PRE(n);
    //                 }
    //             }
    //             if (gomv == RSB_BOOL_TRUE)
    //             {
    //                 const size_t scoff = submatrix->coff - mtxAp->coff;
    //                 const size_t sroff = submatrix->roff - mtxAp->roff;
    //                 const char *const offx = ((const char *)x) + (el_size * scoff) * incx;
    //                 char *const offy = ((char *)ov) + (el_size * sroff) * oincy;

    //                 assert(scoff >= 0);
    //                 assert(sroff >= 0);
    //                 RSB_DO_ERROR_CUMULATE(errval, rsb__do_spmv_non_recursive(submatrix, offx, offy, alphap, NULL, incx, oincy, transA RSB_INNER_NRHS_SPMV_ARGS_IDS));

    //                 // #pragma omp critical(rsb_spmv_crs)
    //                 // {
    //                 //     rsb_do_spmv_lock_release(&lock, th_id, ov);
    //                 //     RSB_DO_SPMV_LOCK_DM_INC(lock);
    //                 // }
    //                 // RSB_SPMV_VS_MARK_POST(n);
    //             }
    //         }
    //         // #pragma omp critical(rsb_spmv_crs)
    //         // {
    //         //     dm = RSB_DO_SPMV_LOCK_DM(lock);
    //         // }

    //         if (dm < all_leaf_matrices_n
    // #if RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV
    //             && ((all_leaf_matrices_n - dm) > th_id)
    // #endif // RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV
    //         )
    //         {
    //             goto again;
    //         }
    //     }
}

void rsb_cuda__do_spmv_recursive(const struct rsb_mtx_t *mtxAp, const void *x, void *y, const void *alphap, const void *betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPMV_ARGS)
{
    rsb_err_t errval = RSB_ERR_NO_ERROR;

    const struct rsb_translated_matrix_t *all_leaf_matrices = NULL;
    const rsb_submatrix_idx_t all_leaf_matrices_n = mtxAp->all_leaf_matrices_n;
    const size_t el_size = mtxAp->el_size;

    // IMPORTANT
    if (!rsb__is_recursive_matrix(mtxAp->flags))
    {
        // dim3 dimBlock(2, 2);
        // dim3 dimGrid(2, 2);
        rsb_cuda__spmv_kernel<<<2, 2>>>(mtxAp, x, y, alphap, betap, incx, incy, transA RSB_INNER_NRHS_SPMV_ARGS_IDS);
        // return rsb__do_spmv_non_recursive(mtxAp, x, y, alphap, betap, incx, incy, transA RSB_INNER_NRHS_SPMV_ARGS_IDS);
    }

    all_leaf_matrices = mtxAp->all_leaf_matrices;

    if (!all_leaf_matrices || all_leaf_matrices_n < 1)
    {
        // errval = RSB_ERR_ENOMEM;
        // goto err;
        return RSB_ERR_ENOMEM;
    }

err:
    return errval;
}

rsb_err_t rsb_cuda__do_spmv(
    rsb_trans_t transA,
    const void *alphap,
    const struct rsb_mtx_t *mtxAp,
    const void *x,
    rsb_coo_idx_t incx,
    const void *betap,
    void *y,
    rsb_coo_idx_t incy,
    enum rsb_op_flags_t op_flags RSB_OUTER_NRHS_SPMV_ARGS)
{

    rsb_err_t errval = RSB_ERR_BADARGS;
    const rsb_type_t typecode = mtxAp->typecode;

#if RSB_ALLOW_ZERO_DIM
    if (RSB_ANY_MTX_DIM_ZERO(mtxAp))
    {
        errval = RSB_ERR_NO_ERROR;
        goto err; /* FIXME: skipping further checks */
    }
#endif

    if (x == y)
    {
        goto err;
    }

    if (incx < 1 || incy < 1)
    {
        goto err;
    }

    if (!alphap || !betap)
    {
        goto err;
    }

    if (!mtxAp || !x || !y || transA == RSB_INVALID_FLAGS)
    {
        goto err;
    }

#if RSB_WANT_OUTER_SPMM_BETA_SCALE
    if (betap && !rsb__is_element_one(betap, typecode))
    {
        RSB_CBLAS_X_SCAL_SPMM(typecode, rsb__do_get_rows_of(mtxAp, transA), betap, y, incy);
        betap = NULL;
    }
#endif /* RSB_WANT_OUTER_SPMM_BETA_SCALE */

#if RSB_WANT_OMP_RECURSIVE_KERNELS
    if (RSB_LIKELY(op_flags != RSB_OP_FLAG_WANT_SERIAL))
    {
        RSB_NUM_THREADS_DECL
        RSB_NUM_THREADS_PUSH
        errval = rsb__do_spmv_recursive_parallel(mtxAp, x, y, alphap, betap, incx, incy, transA, op_flags RSB_OUTER_NRHS_SPMV_ARGS_IDS);
        RSB_NUM_THREADS_POP
    }
    else
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
        errval = rsb_cuda__do_spmv_recursive(mtxAp, x, y, alphap, betap, incx, incy, transA, op_flags RSB_INNER_NRHS_SPMV_ARGS_IDS);

    /* Note: the RSB_OP_FLAG_FAKE_LOCK case is handled by rsb__do_spmv_recursive_parallel */
/*    goto done;
done:
    if (!RSB_UNLIKELY(op_flags & RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT))
    { // NEW: fix for odd spsv/diagonal implicit/no-parallel cases
        if (RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_UNIT_DIAG_IMPLICIT))
        {
            const rsb_coo_idx_t ndy = RSB_MIN(mtxAp->nr, mtxAp->nc);
            const int row_major = (nrhs > 1 && (incx >= nrhs || incy >= nrhs));
            const rsb_nnz_idx_t ldY = rsb__do_get_rows_of(mtxAp, transA);
            const rsb_nnz_idx_t ldX = rsb__do_get_columns_of(mtxAp, transA);
            rsb_int_t nrhsi, di;

            if (row_major)
            {
                for (di = 0; di < ndy; ++di)
                {
                    rsb__cblas_Xaxpy(typecode, nrhs, alphap,
                                     RSB_TYPED_OFF_PTR(typecode, x, incx * di),
                                     1,
                                     RSB_TYPED_OFF_PTR(typecode, y, incy * di), 1);
                }
            }
            else
            {
                for (nrhsi = 0; nrhsi < nrhs; ++nrhsi)
                {
                    rsb__BLAS_Xaxpy_parallel(ndy, alphap, RSB_TYPED_OFF_PTR(typecode, y, incy * nrhsi * ldY), incy, RSB_TYPED_OFF_PTR(typecode, x, incx * nrhsi * ldX), incx, typecode);
                }
            }
        }
    }*/
err:
    RSB_DO_ERR_RETURN(errval)
}

rsb_err_t rsb_cuda_spmv(
    rsb_trans_t transA,
    const void *alphap,
    const struct rsb_mtx_t *mtxAp,
    const void *Xp,
    rsb_coo_idx_t incX,
    const void *betap,
    void *Yp,
    rsb_coo_idx_t incY)
{
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    RSB_INTERFACE_PREAMBLE
    errval = rsb_cuda__do_spmv(transA, alphap, mtxAp, Xp, incX, betap, Yp, incY, (RSB_OP_FLAG_DEFAULT)RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS);
    RSB_INTERFACE_RETURN_ERR(errval)
}
