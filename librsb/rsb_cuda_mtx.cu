#include "rsb_cuda.h"
#include "rsb_do.h"

#include <cuda.h>

#define RSB_SUBDIVISION_BUG_EXTRA (4)

RSB_INTERNALS_COMMON_HEAD_DECLS

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

static rsb_submatrix_idx_t rsb_cuda__estimate_subm_count(const rsb_nnz_idx_t nnz, const rsb_type_t typecode, const rsb_flags_t flags, const rsb_thread_t wet, rsb_err_t *errvalp)
{
	cudaDeviceProp prop;
	int current_device;
	cudaCheckError(cudaGetDevice(&current_device));
	cudaCheckError(cudaGetDeviceProperties(&prop, current_device));

	const long fcs = prop.sharedMemPerBlock;	  // rsb__get_first_level_c_size();
	const long lcs = prop.totalGlobalMem;		  // rsb__get_lastlevel_c_size();
	const long lcspt = prop.totalGlobalMem / nnz; // rsb__get_lastlevel_c_size_per_thread();
	const long cbs = lcspt;						  // rsb__get_cache_block_byte_size();
#if RSB_WANT_SUBDIVISION_FIXES_20101213
	const long wcbs = cbs;
#else							 /* RSB_WANT_SUBDIVISION_FIXES_20101213 */
	const long wcbs = lcspt;
#endif							 /* RSB_WANT_SUBDIVISION_FIXES_20101213 */
	rsb_submatrix_idx_t tmc = 0; /* total (max) matrix count */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if (fcs > lcs || fcs < 1 || cbs < 1) /* we allow declaration of 1 level of cache only */
	{
		/* TODO : find a reasonable solution, and declare it in ./rsbench -I, which should give some diagnostic about this */
		errval = RSB_ERR_FAILED_MEMHIER_DETECTION;
		RSB_PERR_GOTO(err, "innermost cache size:%d, outermost cache size:%d, cache block size %d\n", fcs, lcs, cbs);
	}

	tmc = RSB_SUBDIVISION_BUG_EXTRA + 2 * (((nnz + wcbs) * (rsb__sizeof(typecode) + 2 * sizeof(rsb_coo_idx_t))) / (wcbs)); /* TODO: clean this up */
	tmc = RSB_MAX(1, (rsb_submatrix_idx_t)(rsb_global_session_handle.subdivision_multiplier * tmc));

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS))
		tmc = RSB_MAX(tmc, wet);

#if !RSB_OLD_COO_CRITERIA
	if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_QUAD_PARTITIONING))
	{
		tmc = 1; /* No quad partitioning ? then single submatrix. */
	}
#endif
err:
	RSB_CONDITIONAL_ERRPSET(errvalp, errval);
	RSB_DEBUG_ASSERT(tmc > 0);
	return tmc;
}

struct rsb_mtx_t *rsb_cuda__allocate_recursive_sparse_matrix_from_row_major_coo(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t *pinfop, rsb_flags_t flags, rsb_err_t *errvalp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	struct rsb_mtx_t *submatrices = NULL;
	struct rsb_mtx_t **submatricesp = NULL;
	struct rsb_mtx_t *mtxAp = NULL;
	rsb_time_t hct = RSB_TIME_ZERO;
	rsb_time_t drt = RSB_TIME_ZERO;
	rsb_time_t mat = RSB_TIME_ZERO;
	rsb_time_t lst = RSB_TIME_ZERO;
	rsb_coo_idx_t *IB = NULL;
	rsb_coo_idx_t *IT = NULL;
	rsb_coo_idx_t *IX = NULL;
	rsb_coo_idx_t *WA = NULL;
	rsb_submatrix_idx_t smi = 0;		  /* submatrix index */
	rsb_submatrix_idx_t cmc = 0, omc = 0; /* closed matrices count, open matrices count */
	rsb_submatrix_idx_t lm = 0;			  /* leaf matrices */
	rsb_time_t dt = RSB_TIME_ZERO;
	rsb_time_t eit = RSB_TIME_ZERO;
	rsb_time_t est = RSB_TIME_ZERO;
	rsb_time_t ect = RSB_TIME_ZERO;
	rsb_time_t tat = RSB_TIME_ZERO;
	rsb_time_t sat = RSB_TIME_ZERO;
	printf("Calling rsb_get_num_coo2rec_threads\n");
	// const rsb_thread_t wet = rsb_get_num_coo2rec_threads(); /* want executing threads: */
	const rsb_thread_t wet = nnz; // 1 (global) CUDA thread for each element
	const size_t el_size = rsb__sizeof(typecode);
	rsb_coo_idx_t roff = 0;
	rsb_coo_idx_t coff = 0;
	rsb_nnz_idx_t dnnz = 0;
	printf("Calling rsb__estimate_subm_count\n");
	const rsb_submatrix_idx_t tmc = rsb_cuda__estimate_subm_count(nnz, typecode, flags, wet, &errval); /* total (max) matrix count */

	if (RSB_SOME_ERROR(errval))
	{
		goto err;
	}

	tat = -(dt = rsb_time());
#if RSB_ALLOW_EMPTY_MATRICES
	if (!RSB_HAVE_GOOD_PARMS_FOR_EMPTY(m, k, nnz, flags))
#endif /* RSB_ALLOW_EMPTY_MATRICES */
		if (!RSB_HAVE_GOOD_PARMS_FOR_IN_PLACE_RCSR(m, k, nnz, flags))
		{
			errval = RSB_ERR_BADARGS;
			RSB_PERR_GOTO(err, RSB_ERRM_MDNFARTS);
		}

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_SYMMETRIC) && RSB_DO_FLAG_HAS(flags, RSB_FLAG_HERMITIAN))
	{
		errval = RSB_ERR_BADARGS;
		RSB_PERR_GOTO(err, RSB_ERRM_BFSAH);
	}

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_WANT_COLUMN_MAJOR_ORDER))
	{
		errval = RSB_ERR_UNIMPLEMENTED_YET;
		RSB_PERR_GOTO(err, RSB_ERRM_CMOINIY);
	}

#if 1
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE))
	{
		errval = RSB_ERR_INTERNAL_ERROR;
		RSB_PERR_GOTO(err, "Unexpected RSB_FLAG_FORTRAN_INDICES_INTERFACE flags here!\n");
	}
#else  /* */
	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE)) /* TODO: this is *slow*, speed this up */
		rsb__util_coo_array_from_fortran_indices(IA, nnz, RSB_BOOL_TRUE),
			rsb__util_coo_array_from_fortran_indices(JA, nnz, RSB_BOOL_TRUE),
			RSB_DO_FLAG_DEL(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE);
#endif /* */

#if defined(RSB__VERBOSE_COO2REC) && RSB__VERBOSE_COO2REC != 0
	RSB_STDOUT("Building a matrix with %ld nnz, %ld x %ld\n", (long int)nnz, (long int)m, (long int)k);
#endif

	if (RSB_DO_FLAGS_EXTRACT_STORAGE(flags) == RSB_FLAG_NOFLAGS)
	{
		RSB_DO_FLAG_ADD(flags, RSB_FLAG_DEFAULT_STORAGE_FLAGS);
	}

	/* TODO: may plug *here* upcoming RSB_WANT_FASTER_EXPERIMENTAL_CONSTRUCTOR stuff */

	mat = -dt;
	submatrices = (rsb_mtx_t *)rsb_cuda__calloc(sizeof(struct rsb_mtx_t) * tmc);
	submatricesp = (rsb_mtx_t **)rsb_cuda__calloc(sizeof(struct rsb_mtx_t *) * tmc);

	if (!submatrices || !submatricesp)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}
	mat += (dt = rsb_time());

	for (smi = 0; smi < tmc; ++smi)
		submatricesp[smi] = submatrices + smi;

	ect = -dt;
	if ((errval = rsb__do_cleanup_nnz(VA, IA, JA, nnz, roff, coff, m, k, &nnz, typecode, flags)) != RSB_ERR_NO_ERROR)
		goto err;
	ect += (dt = rsb_time());

	est = -dt;
	if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_SORTED_INPUT))
		if ((errval = rsb__util_sort_row_major_inner(VA, IA, JA, nnz, m, k, typecode, flags)) != RSB_ERR_NO_ERROR)
			RSB_PERR_GOTO(err, RSB_ERRM_ES);

	RSB_DO_FLAG_ADD(flags, RSB_FLAG_SORTED_INPUT); /* TODO: is this needed ? */
	est += (dt = rsb_time());					   /* with 'sorting' (est) we DO NOT intend also cleanup (in ect) */

	/* we need duplicates removal, and this can only take place after sorting */
	drt = -dt;
	dnnz = nnz - rsb__weed_out_duplicates(IA, JA, VA, nnz, typecode, flags);
	nnz -= dnnz;
	drt += (dt = rsb_time());
#if defined(RSB__VERBOSE_COO2REC) && RSB__VERBOSE_COO2REC != 0
	RSB_INFO("Duplicates check: %zd - %zd = %zd\n", (size_t)(nnz + dnnz), (size_t)dnnz, (size_t)nnz);
#endif

	/* work vectors allocation */
	/*	IL = rsb__malloc(sizeof(rsb_coo_idx_t)*(m+1)); */
	mat -= dt;
	IT = (rsb_coo_idx_t *)rsb_cuda__malloc(sizeof(rsb_coo_idx_t) * (m + 1));
	IX = (rsb_coo_idx_t *)rsb_cuda__malloc(sizeof(rsb_coo_idx_t) * 2 * (m + 1));
	IB = (rsb_coo_idx_t *)rsb_cuda__malloc(sizeof(rsb_coo_idx_t) * (m + 1));
	mat += (dt = rsb_time());
	if (/*  !IL ||*/ !IT || !IX || !IB)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}

	/* declaring this first matrix (smi == 0) as 'open' */
	smi = 0;
	omc = 1;
	/* compile in the first mtxAp, linking into to the temporary split vector */
	submatrices[smi].nzoff = 0;
	submatrices[smi].roff = roff;
	submatrices[smi].coff = coff;
	submatrices[smi].bindx = IB;
	submatrices[smi].bpntr = IB + 1;
	submatrices[smi].indptr = NULL;
	/*	RSB_DO_FLAG_ADD(flags,RSB_FLAG_QUAD_PARTITIONING); */
	RSB_DO_FLAG_ADD(flags, RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS);
	if ((errval = rsb__set_init_flags_and_stuff(
			 &submatrices[smi], NULL, NULL, m, k, nnz, nnz, nnz, typecode, flags)) != RSB_ERR_NO_ERROR)
		RSB_PERR_GOTO(err, RSB_ERRM_ES);

	if (nnz == 0)
	{
		/* a special case. we copy the arrays addresses because they may be non-NULL and containing duplicate/diagonal/etc.. elements we have honoured to free, afterwards. */
		++cmc;
		--omc;
		mtxAp = &submatrices[0];
		RSB_DO_FLAG_DEL(mtxAp->flags, RSB_FLAG_QUAD_PARTITIONING); /* necessary, too */
		mtxAp->bpntr = IA;
		mtxAp->bindx = JA;
		mtxAp->indptr = NULL;
		mtxAp->VA = VA;
		goto arrays_done;
	}
	else
		mtxAp = &submatrices[0];

	sat = -(dt = rsb_time());
	/* computing the first right-left pointer vectors */
	if ((errval = rsb_do_compute_vertical_split_parallel(IA, JA, roff, coff, m, k, 0, 0, nnz, IB, NULL, NULL, NULL, NULL, NULL, NULL)) != RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}

	if (RSB_C2R_IF_VERBOSE)
		RSB_INFO("beginning (%zd x %zd) @ %p with flags 0x%x (coo:%d, csr:%d), storage: 0x%x, max %zd submatrices\n",
				 (rsb_printf_int_t)submatrices[smi].nr, (rsb_printf_int_t)submatrices[smi].nc, (const void *)&submatrices[smi], submatrices[smi].flags,
				 RSB_DO_FLAG_HAS(submatrices[smi].flags, RSB_FLAG_WANT_COO_STORAGE),
				 RSB_DO_FLAG_HAS(submatrices[smi].flags, RSB_FLAG_WANT_BCSS_STORAGE),
				 submatrices[smi].matrix_storage, (rsb_printf_int_t)tmc);

/*	if(!RSB_WANT_MORE_PARALLELISM || (RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_QUAD_PARTITIONING))) */ /* TODO */
#if 1																								/* the code is not yet ready for this */
																									/* #if RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = rsb_do_coo2rec_subdivide_parallel(VA, IA, JA, m, k, nnz, typecode, pinfop, flags, errvalp, submatricesp, mtxAp, IB, IX, IT, WA, cmc, omc, tmc, RSB_MAX(1, RSB_MIN(wet, nnz)), &cmc);
#else
	errval = rsb_do_coo2rec_subdivide(VA, IA, JA, m, k, nnz, typecode, pinfop, flags, errvalp, submatricesp, mtxAp, IB, IX, IT, WA, cmc, omc, tmc, wet, &cmc);
#endif
	sat += (dt = rsb_time());

	RSB_CONDITIONAL_FREE(IX);
	if (RSB_SOME_ERROR(errval))
		goto err;

	/*
	RSB_CONDITIONAL_FREE(IL);
	RSB_CONDITIONAL_FREE(IT);
		*/
	/* WA will is needed for shuffle only (so, after IL,IM deallocation, in a way that total memory need is max(WA,IL)) */
	mat -= dt;
	WA = rsb_cuda__malloc(RSB_MAX(sizeof(rsb_coo_idx_t), el_size) * nnz);
	if (!WA)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err, RSB_ERRM_FAOTAFS);
	}
	mat += (dt = rsb_time());

	eit = -dt;

	for (smi = 0; smi < cmc; ++smi)
		if (rsb__is_terminal_recursive_matrix(submatricesp[smi]))
			++lm;

	/*	qsort(submatricesp+(cmc-lm),(size_t)(lm),sizeof(struct rsb_mtx_t*),&rsb__compar_rcsr_matrix_leftmost_first); */
	qsort(submatricesp, (size_t)(cmc), sizeof(struct rsb_mtx_t *), &rsb__compar_rcsr_matrix_leftmost_first);
	/* TODO: a priority queue would do the job, here */
	for (smi = 0; smi < cmc - lm; ++smi)
		if (rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ANLSMIT);
		}
	for (smi = cmc - lm; smi < cmc; ++smi)
		if (!rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ALSMINT);
		}

	errval = rsb_do_coo2rec_shuffle(VA, IA, JA, m, k, nnz, typecode, pinfop, flags, errvalp, submatricesp, mtxAp, IB, WA, cmc);
	if (RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err, RSB_ERRM_SEOWS);
	}

	rsb__do_set_in_place_submatrices_offsets(submatrices, cmc, VA, IA, JA, el_size);

	/*	RSB_INFO("VA:%p, IA:%p, JA:%p\n",VA,IA,JA); */

	if (RSB_C2R_IF_VERBOSE)
		RSB_INFO("IA? :%p / %p\n", (const void *)IA,
				 (const void *)(rsb__do_get_first_submatrix(mtxAp)->bpntr -
								(rsb__do_get_first_submatrix(mtxAp)->nr + 1))
				 /*			rsb__do_get_first_submatrix(mtxAp)->roff-
							 ((submatricesp[0])->nr+1) */
		);

	/* after shuffling, the last vectors conversion happens and we are done. */
arrays_done:
	eit += (dt = rsb_time());

	lst -= dt;
	/* first matrix is always root (even if a CSR one) */
	RSB_DO_FLAG_DEL(submatrices[0].flags, RSB_FLAG_NON_ROOT_MATRIX);
#if RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS
	if (!(submatrices[0].flags & RSB_FLAG_NON_ROOT_MATRIX))
	{
		submatrices[0].all_leaf_matrices = NULL;
		errval = rsb__get_array_of_leaf_matrices(&submatrices[0], &(submatrices[0].all_leaf_matrices), &submatrices[0].all_leaf_matrices_n);
		if (RSB_SOME_ERROR(errval))
			goto err;
	}
	else
	{
		/* this is a non root matrix */
		submatrices[0].all_leaf_matrices = NULL;
		submatrices[0].all_leaf_matrices_n = 0;
	}
#endif /* RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS */
	lst += (dt = rsb_time());

	hct = -dt;
	if (
		/* RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_COO)
			   ||	RSB_DO_FLAG_HAS(flags,RSB_FLAG_USE_HALFWORD_INDICES_CSR)
			   ||	RSB_DO_FLAG_HAS(flags,RSB_FLAG_WANT_COO_STORAGE)*/
		RSB_DO_FLAG_HAS(flags, RSB_FLAG_USE_HALFWORD_INDICES))
#if RSB_WANT_MORE_PARALLELISM
		RSB_DO_ERROR_CUMULATE(errval, rsb_do_switch_fresh_recursive_matrix_to_halfword_storages_parallel(mtxAp));
#else  /* RSB_WANT_MORE_PARALLELISM */
		RSB_DO_ERROR_CUMULATE(errval, rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(mtxAp));
#endif /* RSB_WANT_MORE_PARALLELISM */
	else
	{
		if (RSB_C2R_IF_VERBOSE)
			RSB_INFO("no  RSB_FLAG_USE_HALFWORD_INDICES flag\n");
	}
	RSB_DO_ERROR_CUMULATE(errval, rsb_do_compute_bounded_boxes(mtxAp));

	hct += rsb_time();

	if (RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}
	tat += rsb_time();

	mtxAp->sat = sat;
	mtxAp->ect = ect;
	mtxAp->eit = eit;
	mtxAp->est = est;
	mtxAp->pet = RSB_TIME_ZERO;
	mtxAp->rpt = RSB_TIME_ZERO;
	mtxAp->tat = tat;
	mtxAp->cpt = RSB_TIME_ZERO; /* cpt is contained in sat, so should not be counted here! */

#if defined(RSB__VERBOSE_COO2REC) && RSB__VERBOSE_COO2REC != 0
	RSB_INFO(" converted COO to RSB in %.3le s (%.2lf %%)\n", tat, RSB_PCT(tat, tat));
	RSB_INFO(" analyzed arrays in %.3le s (%.2lf %%)\n", sat, RSB_PCT(sat, tat));
	RSB_INFO(" cleaned-up arrays in %.3le s (%.2lf %%)\n", ect, RSB_PCT(ect, tat));
	RSB_INFO(" deduplicated arrays in %.3le s (%.2lf %%)\n", drt, RSB_PCT(drt, tat));
	RSB_INFO(" sorted arrays in %.3le s (%.2lf %%)\n", est, RSB_PCT(est, tat));
	if (mtxAp->cpt)
		RSB_INFO(" computed partitions in %.3le s (%.2lf %%)\n", mtxAp->cpt, RSB_PCT(mtxAp->cpt, tat));
	RSB_INFO(" shuffled partitions in %.3le s (%.2lf %%)\n", eit, RSB_PCT(eit, tat));
	RSB_INFO(" memory allocations took %.3le s (%.2lf %%)\n", mat, RSB_PCT(mat, tat));
	RSB_INFO(" leafs setup took %.3le s (%.2lf %%)\n", lst, RSB_PCT(lst, tat));
	RSB_INFO(" halfword conversion took %.3le s (%.2lf %%)\n", hct, RSB_PCT(hct, tat));
#endif

#if RSB_STORE_IDXSA
	mtxAp->idxsa = rsb__get_index_storage_amount(mtxAp);
#endif

	goto noerr;
err:
	mtxAp = NULL;
	RSB_CONDITIONAL_FREE(submatrices);
noerr:
	if (RSB_SOME_ERROR(errval))
		rsb__do_perror(NULL, errval);
	RSB_CONDITIONAL_FREE(IB);
	RSB_CONDITIONAL_FREE(IT);
	RSB_CONDITIONAL_FREE(IX);
	/*	RSB_CONDITIONAL_FREE(IL); */
	RSB_CONDITIONAL_FREE(WA);
	RSB_CONDITIONAL_FREE(submatricesp);
	RSB_CONDITIONAL_ERRPSET(errvalp, errval);

#if defined(RSB__VERBOSE_COO2REC) && RSB__VERBOSE_COO2REC != 0
	if (mtxAp)
	{
		RSB_STDOUT("Built ");
		RSB_STDOUT_MATRIX_SUMMARY(mtxAp);
		RSB_STDOUT("\n");
	}
#endif
	return mtxAp;
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
			printf("Calling rsb__util_find_coo_max_index_val\n");
			m = rsb__util_find_coo_max_index_val(IA, nnz) + roff + 1;
		}
		if (k == 0 && JA)
		{
			printf("Calling rsb__util_find_coo_max_index_val\n");
			k = rsb__util_find_coo_max_index_val(JA, nnz) + coff + 1;
		}
		// printf("rc %d %d %d \n",m,k,nnz);
	}

	if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_FORTRAN_INDICES_INTERFACE))
	{
		RSB_PERR_GOTO(err, "!\n");
	}

	if (RSB__FLAG_HAS_UNSPECIFIED_TRIANGLE(flags))
	{
		printf("Calling rsb__do_detect_and_add_triangular_flags\n");
		RSB_DO_FLAG_ADD(flags, rsb__do_detect_and_add_triangular_flags(IA, JA, nnz, flags));
	}

	if (roff && IA)
	{
		printf("Calling rsb__util_coo_array_add\n");
		rsb__util_coo_array_add(IA, nnz, roff);
	}
	if (coff && JA)
	{
		printf("Calling rsb__util_coo_array_add\n");
		rsb__util_coo_array_add(JA, nnz, coff);
	}

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
			printf("Calling rsb__allocate_recursive_sparse_matrix_from_row_major_coo\n");
			mtxAp = rsb_cuda__allocate_recursive_sparse_matrix_from_row_major_coo(VA, IA, JA, m, k, nnz, typecode, NULL, flags, errvalp);
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
		printf("Calling rsb_cuda__util_coo_alloc_copy_and_stats\n");
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

	printf("Calling rsb_cuda__mtx_alloc_inner\n");
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
	printf("calling rsb_cuda__do_mtx_alloc_from_coo_const\n");
	mtxAp = rsb_cuda__do_mtx_alloc_from_coo_const(VA, IA, JA, nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp, errval, errvalp);
}

struct rsb_mtx_t *rsb_cuda_mtx_free(struct rsb_mtx_t *mtxAp)
{
	printf("Called rsb_cuda_mtx_free which is not implemented\n");
	return NULL;
}