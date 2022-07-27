#include "rsb_cuda.h"
#include "rsb_do.h"

#define RSB_INTERFACE_RETURN_MTX_ERRP(MTXAP, ERRVAL, ERRVALP) \
    RSB_INTERFACE_ENDCMD                                      \
    RSB_CONDITIONAL_ERRPSET(ERRVALP, ERRVAL)                  \
    RSB_DO_MTX_RETURN_INTERFACE(MTXAP, ERRVAL);

#define RSB_INITIALIZE_CHECK_MTX_ERRP(ERRVALP)                                      \
    if (!rsb__do_was_initialized())                                                 \
    {                                                                               \
        RSB_ERROR(RSB_ERRM_UL);                                                     \
        RSB_INTERFACE_RETURN_MTX_ERRP(NULL, RSB_ERR_UNSUPPORTED_OPERATION, ERRVALP) \
    }

struct rsb_mtx_t * rsb_cuda__do_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flags, rsb_err_t * errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t * mtxAp = NULL;
	blas_sparse_matrix bmtxA = RSB_BLAS_INVALID_VAL;

	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_SYMMETRIC) && RSB_DO_FLAG_HAS(flags,RSB_FLAG_HERMITIAN))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_BFSAH);
	}

    cudaCheckError( cudaMallocManaged(&mtxAp, sizeof(struct rsb_mtx_t)) );
    rsb__init_struct(mtxAp);
	if(!mtxAp)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAP"\n");
	}
	RSB_MTX_SET_HBDF(mtxAp);
	bmtxA = mtxAp->RSB_MTX_BDF = rsb__BLAS_Xuscr_begin(nrA,ncA,typecode);
	if( mtxAp->RSB_MTX_BDF == RSB_BLAS_INVALID_VAL )
	{
		errval = RSB_ERR_GENERIC_ERROR;
		RSB_CONDITIONAL_FREE(mtxAp);
		RSB_PERR_GOTO(err,RSB_ERRM_IPEWIEM);
	}

	/* FIXME : the following need an improvement  */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_FORTRAN_INDICES_INTERFACE)) rsb__BLAS_ussp( bmtxA, blas_one_base);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UNIT_DIAG_IMPLICIT)) rsb__BLAS_ussp( bmtxA, blas_unit_diag );
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_TRIANGULAR)) rsb__BLAS_ussp( bmtxA, blas_lower_triangular);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_TRIANGULAR)) rsb__BLAS_ussp( bmtxA, blas_upper_triangular);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_TRIANGULAR)) rsb__BLAS_ussp( bmtxA, blas_triangular); /* ask for detection */
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_SYMMETRIC)) rsb__BLAS_ussp( bmtxA, blas_lower_symmetric);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_SYMMETRIC)) rsb__BLAS_ussp( bmtxA, blas_upper_symmetric);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_LOWER_HERMITIAN)) rsb__BLAS_ussp( bmtxA, blas_lower_hermitian);
	if(RSB_DO_FLAG_HAS(flags,RSB_FLAG_UPPER_HERMITIAN)) rsb__BLAS_ussp( bmtxA, blas_upper_hermitian);
err:
	RSB_CONDITIONAL_ERRPSET(errvalp,errval);
	return mtxAp;
}

rsb_err_t rsb_cuda__do_mtx_alloc_from_coo_end(struct rsb_mtx_t ** mtxApp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	blas_sparse_matrix bmtxA = RSB_BLAS_INVALID_VAL;
	struct rsb_mtx_t * mtxBp = NULL;
	struct rsb_mtx_t * mtxAp = NULL;

	if(!mtxApp || !*mtxApp)
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_E_MTXAPP);
	}

	mtxAp = *mtxApp ;

	if( !RSB_MTX_HBDF( mtxAp ) )
	{
		/* errval = RSB_ERR_NO_ERROR; */
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_DNSAMIWAFCB);
	}

	bmtxA = RSB_MTX_HBDFH(mtxAp);
	/* FIXME: missing serious check on mtxAp->flags ! */
	if( rsb__BLAS_Xuscr_end_flagged(bmtxA,NULL) == RSB_BLAS_INVALID_VAL )
	{
	       	errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err,RSB_ERRM_PFTM);
		/* FIXME: insufficient cleanup */
	}
	mtxBp = rsb__BLAS_inner_matrix_retrieve(bmtxA);
	*mtxApp = mtxBp;
	//rsb__free(mtxAp);
    cudaCheckError( cudaFree(mtxAp) );
	rsb__BLAS_handle_free(bmtxA); /* ignoring return value ... */
err:
	return errval;
}

struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_coo_begin(
    rsb_nnz_idx_t nnzA,
    rsb_type_t typecode,
    rsb_coo_idx_t nrA,
    rsb_coo_idx_t ncA,
    rsb_flags_t flagsA,
    rsb_err_t *errvalp)
{
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    struct rsb_mtx_t *mtxAp = NULL;
    RSB_INTERFACE_PREAMBLE
    RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
    mtxAp = rsb_cuda__do_mtx_alloc_from_coo_begin(nnzA, typecode, nrA, ncA, flagsA, &errval);
    RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp, errval, errvalp);
}

rsb_err_t rsb_cuda_mtx_alloc_from_coo_end(struct rsb_mtx_t **mtxApp)
{
    rsb_err_t errval = RSB_ERR_BADARGS;
    RSB_INTERFACE_PREAMBLE
    errval = rsb_cuda__do_mtx_alloc_from_coo_end(mtxApp);
    RSB_INTERFACE_RETURN_ERR(errval)
}