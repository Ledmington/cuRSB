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

#define RSB_CUDA_CONDITIONAL_FREE(p)       \
	{                                      \
		if ((p))                           \
		{                                  \
			cudaCheckError(cudaFree((p))); \
			(p) = NULL;                    \
		}                                  \
	}

static inline size_t rsb__sizeof(rsb_type_t type)
{
	size_t so = 0;
	switch (type)
	{
	case RSB_NUMERICAL_TYPE_DOUBLE:
		so = sizeof(double);
		break;
	case RSB_NUMERICAL_TYPE_FLOAT:
		so = sizeof(float);
		break;
	case RSB_NUMERICAL_TYPE_FLOAT_COMPLEX:
		so = sizeof(float _Complex);
		break;
	case RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX:
		so = sizeof(double _Complex);
		break;
	/* unsupported type */
	default:
		so = 0;
		break;
	}
	return so;
}

struct rsb_mtx_t *rsb_cuda__do_mtx_alloc_from_coo_begin(rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_flags_t flags, rsb_err_t *errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;
	blas_sparse_matrix bmtxA = RSB_BLAS_INVALID_VAL;

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_SYMMETRIC) && RSB_DO_FLAG_HAS(flags, RSB_FLAG_HERMITIAN))
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	cudaCheckError(cudaMallocManaged(&mtxAp, sizeof(struct rsb_mtx_t)));
	rsb__init_struct(mtxAp);
	if (!mtxAp)
	{
		errval = RSB_ERR_ENOMEM;
		goto err;
	}
	RSB_MTX_SET_HBDF(mtxAp);
	bmtxA = mtxAp->RSB_MTX_BDF = rsb__BLAS_Xuscr_begin(nrA, ncA, typecode);
	if (mtxAp->RSB_MTX_BDF == RSB_BLAS_INVALID_VAL)
	{
		errval = RSB_ERR_GENERIC_ERROR;
		RSB_CUDA_CONDITIONAL_FREE(mtxAp);
		goto err;
	}

	/* FIXME : the following need an improvement  */
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE))
		rsb__BLAS_ussp(bmtxA, blas_one_base);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_UNIT_DIAG_IMPLICIT))
		rsb__BLAS_ussp(bmtxA, blas_unit_diag);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_LOWER_TRIANGULAR))
		rsb__BLAS_ussp(bmtxA, blas_lower_triangular);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_UPPER_TRIANGULAR))
		rsb__BLAS_ussp(bmtxA, blas_upper_triangular);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_TRIANGULAR))
		rsb__BLAS_ussp(bmtxA, blas_triangular); /* ask for detection */
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_LOWER_SYMMETRIC))
		rsb__BLAS_ussp(bmtxA, blas_lower_symmetric);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_UPPER_SYMMETRIC))
		rsb__BLAS_ussp(bmtxA, blas_upper_symmetric);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_LOWER_HERMITIAN))
		rsb__BLAS_ussp(bmtxA, blas_lower_hermitian);
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_UPPER_HERMITIAN))
		rsb__BLAS_ussp(bmtxA, blas_upper_hermitian);
err:
	RSB_CONDITIONAL_ERRPSET(errvalp, errval);
	return mtxAp;
}

rsb_err_t rsb_cuda__do_mtx_alloc_from_coo_end(struct rsb_mtx_t **mtxApp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	blas_sparse_matrix bmtxA = RSB_BLAS_INVALID_VAL;
	struct rsb_mtx_t *mtxBp = NULL;
	struct rsb_mtx_t *mtxAp = NULL;

	if (!mtxApp || !*mtxApp)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	mtxAp = *mtxApp;

	if (!RSB_MTX_HBDF(mtxAp))
	{
		/* errval = RSB_ERR_NO_ERROR; */
		errval = RSB_ERR_BADARGS;
		goto err;
	}

	bmtxA = RSB_MTX_HBDFH(mtxAp);
	/* FIXME: missing serious check on mtxAp->flags ! */
	if (rsb__BLAS_Xuscr_end_flagged(bmtxA, NULL) == RSB_BLAS_INVALID_VAL)
	{
		errval = RSB_ERR_BADARGS;
		goto err;
		/* FIXME: insufficient cleanup */
	}
	mtxBp = rsb__BLAS_inner_matrix_retrieve(bmtxA);
	*mtxApp = mtxBp;
	// rsb__free(mtxAp);
	cudaCheckError(cudaFree(mtxAp));
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

void *rsb_cuda__malloc(size_t size)
{
	void *ptr = NULL;
	cudaCheckError(cudaMallocManaged(&ptr, size));
	return ptr;
}

void *rsb_cuda__calloc(size_t size)
{
	void *ptr = rsb__malloc(size);
	if (ptr)
	{
		RSB_BZERO(ptr, size);
	}
	return ptr;
}

void *rsb_cuda__calloc_vector(rsb_nnz_idx_t n, rsb_type_t typecode)
{
	const size_t so = rsb__sizeof(typecode);
	return rsb_cuda__calloc(so * n);
}

void *rsb_cuda__malloc_vector(rsb_nnz_idx_t n, rsb_type_t typecode)
{
	const size_t so = rsb__sizeof(typecode);
	return rsb_cuda__malloc(so * n);
}

rsb_err_t rsb_cuda__util_coo_alloc(void **RSB_RESTRICT VAp, rsb_coo_idx_t **RSB_RESTRICT IAp, rsb_coo_idx_t **RSB_RESTRICT JAp, rsb_nnz_idx_t nnz, rsb_type_t typecode, rsb_bool_t do_calloc)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL, *IA_ = NULL, *JA_ = NULL;

	if (RSB_MATRIX_UNSUPPORTED_TYPE(typecode))
	{
		errval = RSB_ERR_UNSUPPORTED_TYPE;
		goto err;
	}

	if (do_calloc != RSB_BOOL_TRUE)
	{
		VA_ = rsb_cuda__malloc_vector((nnz), typecode),
		IA_ = rsb_cuda__malloc(sizeof(rsb_coo_idx_t) * (nnz)),
		JA_ = rsb_cuda__malloc(sizeof(rsb_coo_idx_t) * (nnz));
	}
	else
	{
		VA_ = rsb_cuda__calloc_vector((nnz), typecode),
		IA_ = rsb_cuda__calloc(sizeof(rsb_coo_idx_t) * (nnz)),
		JA_ = rsb_cuda__calloc(sizeof(rsb_coo_idx_t) * (nnz));
	}

	if (!VA_ || !IA_ || !JA_)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err, RSB_ERRM_ENOMEM);
	}

	*VAp = VA_;
	*IAp = (rsb_coo_idx_t *)IA_;
	*JAp = (rsb_coo_idx_t *)JA_;
	goto done;
err:
	RSB_CUDA_CONDITIONAL_FREE(IA_);
	RSB_CUDA_CONDITIONAL_FREE(JA_);
	RSB_CUDA_CONDITIONAL_FREE(VA_);
done:
	return errval;
}

rsb_err_t rsb_cuda__util_coo_alloc_copy_and_stats(void **RSB_RESTRICT VAp, rsb_coo_idx_t **RSB_RESTRICT IAp, rsb_coo_idx_t **RSB_RESTRICT JAp, const void *RSB_RESTRICT VA, const rsb_coo_idx_t *RSB_RESTRICT IA, const rsb_coo_idx_t *RSB_RESTRICT JA, rsb_coo_idx_t *RSB_RESTRICT mp, rsb_coo_idx_t *RSB_RESTRICT kp, rsb_nnz_idx_t nnz, rsb_nnz_idx_t ennz, rsb_type_t typecode, const rsb_coo_idx_t offi, const rsb_coo_idx_t offo, rsb_flags_t iflags, rsb_flags_t *RSB_RESTRICT flagsp)
{
	/*!
	 * Copies contents of a COO arrays triple to a freshly allocated COO arrays triple.
	 * Size is assumed to be nnz+ennz.
	 * Last ennz elements are not zeroed.
	 *
	 * Flags are determined: RSB_FLAG_UPPER_TRIANGULAR, RSB_FLAG_LOWER_TRIANGULAR.
	 *
	 * TODO: May implement input sanitization or zeroes detection.
	 * TODO: Check for nnz+ennz overflow.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL;
	rsb_coo_idx_t *IA_ = NULL, *JA_ = NULL;

	errval = rsb_cuda__util_coo_alloc((void **)(&VA_), &IA_, &JA_, nnz + ennz, typecode, RSB_BOOL_FALSE);
	if (RSB_SOME_ERROR(errval))
		goto err;

	if (!VA && !IA && !JA)
		goto nocopy; /* it's ok: alloc only semantics */
	/* TODO: the following shall be made parallel */
	if (mp || kp || (flagsp && RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(*flagsp)))
		errval = rsb__util_coo_copy_and_stats(VA, IA, JA, VA_, IA_, JA_, mp, kp, nnz, typecode, offi, offo, iflags, flagsp);
	else
	{
		errval = rsb__util_coo_copy(VA, IA, JA, VA_, IA_, JA_, nnz, typecode, offi, offo);
		/* ... flags may not always be desired! */
		/*	if(flagsp)
				(*flagsp)|=rsb__util_coo_determine_uplo_flags(IA_,JA_,nnz);*/
	}
nocopy:
	*VAp = VA_;
	*IAp = IA_;
	*JAp = JA_;
	goto done;
err:
	RSB_CUDA_CONDITIONAL_FREE(IA_);
	RSB_CUDA_CONDITIONAL_FREE(JA_);
	RSB_CUDA_CONDITIONAL_FREE(VA_);
done:
	return errval;
}

struct rsb_mtx_t *rsb_cuda__mtx_alloc_inner(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_nnz_idx_t nnz, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_type_t typecode, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_blk_idx_t Mb, rsb_blk_idx_t Kb, rsb_flags_t flags, rsb_err_t *errvalp)
{
	/*!

	   Allocates a blocked sparse matrix being recursively partitioned,
	   but in a data structure which is specified by flags, and
	   thus not necessarily with exact BCSR/BCSC leaves.

	   \return a valid matrix pointer or NULL
	*/

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;

	RSB_DEBUG_ASSERT(roff >= -1 && coff >= -1); /* for Fortran */

	if ((m == 0 || k == 0) && nnz > 0)
	{
		/* as a special case, we detect the m and k boundaries, if nnz>0 and m or k are zero */
		/* TODO: shall use rsb__util_coo_alloc_copy_and_stats instead */
		if (m == 0 && IA)
		{
			m = rsb__util_find_coo_max_index_val(IA, nnz) + roff + 1;
		}
		if (k == 0 && JA)
		{
			k = rsb__util_find_coo_max_index_val(JA, nnz) + coff + 1;
		}
		// printf("rc %d %d %d \n",m,k,nnz);
	}

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE))
	{
		RSB_PERR_GOTO(err, "!\n");
	}

	if (RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(flags))
		RSB_DO_FLAG_ADD(flags, rsb__do_detect_and_add_triangular_flags(IA, JA, nnz, flags));

	if (roff && IA)
		rsb__util_coo_array_add(IA, nnz, roff);
	if (coff && JA)
		rsb__util_coo_array_add(JA, nnz, coff);

	RSB_DO_FLAG_ADD(flags, RSB_FLAG_SORT_INPUT);
	RSB_DO_FLAG_ADD(flags, RSB_FLAG_OWN_PARTITIONING_ARRAYS); /* this is in order to free p_r and p_c with the matrix itself, and ignore original flag on this topic */
	if (
		(m < RSB_MIN_MATRIX_DIM || k < RSB_MIN_MATRIX_DIM || nnz < RSB_MIN_MATRIX_NNZ) ||
		(m > RSB_MAX_MATRIX_DIM || k > RSB_MAX_MATRIX_DIM || nnz > RSB_MAX_MATRIX_NNZ))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}
#if RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT
	if (RSB_DO_TOOFEWNNZFORRCSR(nnz, RSB_MIN(m, k)))
		RSB_DO_FLAG_ADD(flags, RSB_FLAG_WANT_COO_STORAGE);
#endif /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */
	if (RSB_HAVE_GOOD_PARMS_FOR_IN_PLACE_RCSR(m, k, nnz, flags)
#if RSB_ALLOW_EMPTY_MATRICES
		|| (nnz == 0)
#endif /* RSB_ALLOW_EMPTY_MATRICES */
	)
	{
		if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_NON_ROOT_MATRIX))
		{
			mtxAp = rsb__allocate_recursive_sparse_matrix_from_row_major_coo(VA, IA, JA, m, k, nnz, typecode, NULL, flags, errvalp);
			if (errvalp && RSB_SOME_ERROR(*errvalp))
				RSB_ERROR("%s\n", rsb__get_errstr_ptr(*errvalp));
			return mtxAp;
		}
	}
#if RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT
	errval = RSB_ERR_INTERNAL_ERROR;
	RSB_ERROR("Trying to call unsupported parameters combination (nr:%ld nc:%ld nnz:%ld)!\n", (long int)m, (long int)k, (long int)nnz);
	rsb__debug_print_flags(flags);
	RSB_PERR_GOTO(err, RSB_ERRM_INTERNAL_ERROR);
#endif /* RSB_WANT_RSB_AS_ONLY_ALLOWED_FORMAT */

	if (mtxAp)
		return mtxAp;
	else
		goto err;
err:
	RSB_CONDITIONAL_ERRPSET(errvalp, errval);
	return NULL;
}

struct rsb_mtx_t *rsb_cuda__do_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flags, rsb_err_t *errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	void *VA_ = NULL;
	rsb_coo_idx_t *IA_ = NULL, *JA_ = NULL;
	struct rsb_mtx_t *mtxAp = NULL;

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
	{
		errval = /*RSB_ERR_BADARGS|*/ RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS;
		RSB_PERR_GOTO(err, RSB_ERRM_CNHEAF);
	}

	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);

	if (nnzA > 0)
	{
		rsb_coo_idx_t offi = 0;

		if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE))
			offi = 1, RSB_DO_FLAG_DEL(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE);
		errval = rsb_cuda__util_coo_alloc_copy_and_stats(&VA_, &IA_, &JA_, VA, IA, JA, nrA ? NULL : &nrA, ncA ? NULL : &ncA, nnzA, 0, typecode, offi, 0, RSB_FLAG_NOFLAGS, &flags);

		if (!VA_ || !IA_ || !JA_)
		{
			errval = RSB_ERR_ENOMEM;
			RSB_PERR_GOTO(err, RSB_ERRM_E_VIJ);
		}
	}
	else
	{
#if !RSB_ALLOW_EMPTY_MATRICES
		/* FIXUP CASE FOR 0-NNZ MATRICES AND IMPLICIT DIAGONAL */
		if (RSB_INVALID_NNZ_COUNT_FOR_FLAGS(nnzA, flags))
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err, RSB_ERRM_CBAEM);
		}
#endif /* RSB_ALLOW_EMPTY_MATRICES */
	}
	RSB_IF_NOFLAGS_SET_DEFAULT_MATRIX_FLAGS(flags);

	mtxAp = rsb_cuda__mtx_alloc_inner(VA_, IA_, JA_, nnzA, 0, 0, typecode, nrA, ncA, brA, bcA, flags, &errval);
	if (mtxAp && errval == RSB_ERR_NO_ERROR)
		goto ok;
	/* FIXME: and if !matrix but errval ? */
err:
	RSB_CUDA_CONDITIONAL_FREE(IA_);
	RSB_CUDA_CONDITIONAL_FREE(JA_);
	RSB_CUDA_CONDITIONAL_FREE(VA_);
ok:
	RSB_CONDITIONAL_ERRPSET(errvalp, errval);
	return mtxAp;
}

struct rsb_mtx_t *rsb_cuda_mtx_alloc_from_coo_const(const void *VA, const rsb_coo_idx_t *IA, const rsb_coo_idx_t *JA, rsb_nnz_idx_t nnzA, rsb_type_t typecode, rsb_coo_idx_t nrA, rsb_coo_idx_t ncA, rsb_blk_idx_t brA, rsb_blk_idx_t bcA, rsb_flags_t flagsA, rsb_err_t *errvalp)
{

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *mtxAp = NULL;
	RSB_INTERFACE_PREAMBLE
	RSB_INITIALIZE_CHECK_MTX_ERRP(errvalp);
	mtxAp = rsb_cuda__do_mtx_alloc_from_coo_const(VA, IA, JA, nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp, errval, errvalp);
}

struct rsb_mtx_t *rsb_cuda_mtx_free(struct rsb_mtx_t *mtxAp)
{
	printf("Called rsb_cuda_mtx_free which is not implemented\n");
	return NULL;
}