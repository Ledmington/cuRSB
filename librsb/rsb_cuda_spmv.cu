#include "rsb_cuda.h"
#include "rsb_internals.h"
#include "rsb_do.h"
//#include "rsb_common.h"

#define RSB_INTERFACE_RETURN_ERR(ERRVAL) RSB_INTERFACE_ENDCMD RSB_DO_ERR_RETURN_INTERFACE(ERRVAL)

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
        goto err;

    if (incx < 1 || incy < 1)
        goto err;

    if (!alphap || !betap)
        goto err;

    if (!mtxAp || !x || !y || transA == RSB_INVALID_FLAGS)
        goto err;

#if RSB_WANT_OUTER_SPMM_BETA_SCALE
    if (betap && !RSB_IS_ELEMENT_ONE(betap, typecode))
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
        errval = rsb__do_spmv_recursive_serial(mtxAp, x, y, alphap, betap, incx, incy, transA RSB_INNER_NRHS_SPMV_ARGS_IDS);

    /* Note: the RSB_OP_FLAG_FAKE_LOCK case is handled by rsb__do_spmv_recursive_parallel */
    goto done;
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
                    1, RSB_TYPED_OFF_PTR(typecode, y, incy * di), 1);
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
    }
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
