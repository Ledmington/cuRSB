#include "rsb_cuda.h"
#include "rsb_do.h"
#include "rsb_common.h"
#include "rsb_coo2rec.h"

#include <cuda.h>

#define RSB_SUBDIVISION_BUG_EXTRA (4)
#define RSB_WANT_BINSEARCH_MIN_NZPR 8
#define RSB_WANT_ZERO_ON_DESTROY 0

#define RSB_A_MEMCPY_SMALL(ID, IS, DOFF, SOFF, NNZ, ES) RSB_A_MEMCPY(ID, IS, DOFF, SOFF, NNZ, ES)
#define RSB_COA_MEMCPY_SMALL(ID, IS, DOFF, SOFF, NNZ) RSB_COA_MEMCPY(ID, IS, DOFF, SOFF, NNZ)

#define RSB_DO_ENOUGHNNZFORINDEXBASEDBUILD(M) (!RSB_DO_TOOFEWNNZFORRCSR((M)->nnz, (M)->nr))

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
		if ((p) != NULL)                   \
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
	void *ptr = rsb_cuda__malloc(size);
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

static rsb_err_t rsb_cuda_do_compute_vertical_split_parallel(const rsb_coo_idx_t *RSB_RESTRICT IA, const rsb_coo_idx_t *RSB_RESTRICT JA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_coo_idx_t hm, rsb_coo_idx_t hk, rsb_nnz_idx_t nnz, rsb_coo_idx_t *IL, rsb_coo_idx_t *RSB_RESTRICT IM, rsb_coo_idx_t *IR, rsb_nnz_idx_t *ulp, rsb_nnz_idx_t *urp, rsb_nnz_idx_t *llp, rsb_nnz_idx_t *lrp)
{
	/**
	Binary search for the boundaries of each row.
	Assign threads to rows intervals.
	Perform the row pointers vector fill calling rsb_do_compute_vertical_split.
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	/* const rsb_thread_t wet = rsb_get_num_threads(); */

	if (m < 1)
	{
		return RSB_ERR_NO_ERROR; /* TODO: limit case */
	}
	IL[0] = 0;
	if (m == 1)
	{
		goto after;
	}

#pragma omp parallel RSB_NTC
	{
#if RSB_WANT_OMP_RECURSIVE_KERNELS
		const rsb_thread_t tn = /*wet*/ omp_get_num_threads(), tnn = RSB_MIN(tn, m);
		const rsb_thread_t th_id = omp_get_thread_num();
#else  /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		const rsb_thread_t tn = 1, tnn = RSB_MIN(tn, m);
		const rsb_thread_t th_id = 0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
		const rsb_coo_idx_t mm = ((m + tnn - 1) / tnn), m0 = mm * th_id, m1 = RSB_MIN(m0 + mm, m);
		rsb_coo_idx_t i;
		rsb_nnz_idx_t nnz0 = 0, nnz1 = nnz;
		rsb_nnz_idx_t n, rnnz, fnnz, lnnz;

		if (th_id >= m)
			goto nowork;
		/* binary search for the boundaries of each row  */
		nnz0 = rsb__nnz_split_coo_bsearch(IA + nnz0, m0, nnz1 - nnz0);
		nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(IA + nnz0, m1, nnz1 - nnz0);
		/* assign threads to rows intervals */
		if (nnz0 >= nnz1)
		{
			// for(i=m0;RSB_LIKELY(i<m1);++i)
			//	IL[i+1] = nnz0;
			RSB_XCOO_VSET(IL, nnz0, m0 + 1, m1 + 1);
			goto nowork;
		}
		// RSB_INFO("thread %d  rows %d..%d  nnz %d..%d\n",th_id,m0,m1,nnz0,nnz1);
		/* perform the row pointers vector fill calling rsb_do_compute_vertical_split */
		// RSB_DO_ERROR_CUMULATE(errval,rsb_do_compute_vertical_split(IA+nnz0,JA+nnz0,roff+m0,coff,m1-m0,k,hm,hk,nnz1-nnz0,IL+m0,NULL,NULL,ulp,urp,llp,lrp));
		fnnz = nnz0;
		n = nnz0;
		for (i = m0; RSB_LIKELY(i < m1); ++i)
		{
			if ((nnz1 - nnz0) / (m1 - m0) < RSB_WANT_BINSEARCH_MIN_NZPR)
			{
				rnnz = 0;
				for (; RSB_LIKELY(n < nnz1 && IA[n] == i); ++n)
					++rnnz;
				IL[i + 1] = nnz0 + rnnz;
				nnz0 += rnnz;
			}
			else
			{
				/* TODO : should use a smarter strategy than this one */
				lnnz = fnnz + rsb__nnz_split_coo_bsearch(IA + fnnz, i + 1, nnz1 - fnnz);
				// RSB_INFO("%d : %d\n",i,lnnz);
				IL[i + 1] = lnnz;
				fnnz = lnnz;
			}
		}
	nowork:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
#pragma omp barrier
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
	}
after:
	IL[m] = nnz;
	// int i; RSB_INFO(":::"); for(i=0;RSB_LIKELY(i<m+1);++i) RSB_INFO("%d ",IL[i]); RSB_INFO("\n");
	// RSB_INFO(":::"); for(i=0;RSB_LIKELY(i<nnz);++i) RSB_INFO("%d ",IA[i]); RSB_INFO("\n");
	// err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_submatrix_idx_t rsb_do_pick_largest_open_matrix(struct rsb_mtx_t **submatricesp, rsb_submatrix_idx_t smc)
{
	rsb_submatrix_idx_t smi = 0, msmi = RSB_SUBM_IDX_MARKER;
	rsb_nnz_idx_t maxnz = 0;
#if defined(RSB_WANT_VERBOSE_SUBDIVISION) && RSB_WANT_VERBOSE_SUBDIVISION != 0
	if (smc == 0)
		RSB_INFO("warning: no largest open matrix among 0 matrices\n");
#endif
	for (smi = 0; smi < smc; ++smi)
	{
		// RSB_INFO("looking %d : %d\n",smi,submatricesp[smi]->nnz);
		/* NOTE: ">=" down here is used to cope with diagonal implicit matrices (which could have nnz==0), too */
		if (submatricesp[smi])
			if (submatricesp[smi]->nnz >= maxnz)
			{
				maxnz = submatricesp[smi]->nnz;
				msmi = smi;
			}
	}
	return msmi;
}

static rsb_bool_t rsb__should_recursively_partition_matrix(
	rsb_coo_idx_t mB, rsb_coo_idx_t kB,
	rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_nnz_idx_t element_count,
	rsb_nnz_idx_t block_count,
	rsb_nnz_idx_t nnz,
	rsb_blk_idx_t Mdim,
	rsb_blk_idx_t mdim,
	rsb_coo_idx_t roff,
	rsb_coo_idx_t coff,
	rsb_flags_t flags,
	size_t el_size,
	rsb_thread_t wet)
{
#if (RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY == RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE)
	long cs = (rsb__get_lastlevel_c_size() / (wet > 0 ? wet : rsb_get_num_threads()));
	rsb_bool_t sp = RSB_BOOL_FALSE; /* should partition */
	rsb_fillin_t efillin = 1.0;		/* FIXME */
	size_t smab = 0;				/* spmv memory accessed bytes */

	if (nnz < RSB_RECURSION_MIN_NNZ || m < RSB_RECURSION_MIN_DIM || k < RSB_RECURSION_MIN_DIM || !RSB_DO_FLAG_HAS(flags, RSB_FLAG_QUAD_PARTITIONING))
	{
		sp = RSB_BOOL_FALSE;
		goto done;
	}

	if (kB < 1)
		kB = 1;
	if (mB < 1)
		mB = 1;

	if ((flags & RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG) && roff == coff)
		cs /= 2;
	/* this will imply a more fine grained subdivision on the diagonal
	 * FIXME : we could use a factor different than 2 !
	 * */

	/* subdivide at least until matrix indices can be compressed */
	if ((flags & RSB_FLAG_USE_HALFWORD_INDICES_CSR) && m > 1 && k > 1 &&
		//			nnz>1
		//			nnz>(cs/4)
		nnz * el_size > 2 * cs && !rsb__do_is_candidate_size_for_halfword_csr(m, k, nnz, flags))
		return RSB_BOOL_TRUE;
	if ((flags & RSB_FLAG_USE_HALFWORD_INDICES_COO) && m > 1 && k > 1 &&
		//			nnz>1
		//			nnz>(cs/4)
		nnz * el_size > 2 * cs && !rsb__do_is_candidate_size_for_halfword_coo(m, k, flags))
		return RSB_BOOL_TRUE;

	if (cs > 0)
	{
		smab = rsb_spmv_memory_accessed_bytes_(mB, kB, m, k, efillin * nnz, ((efillin * nnz) / mB) / kB, m / mB, el_size);

		if (2 * smab > 3 * 3 * cs) /* FIXME : overflow possible */
			sp = 1;
		else if (
			/* FIXME! */
			(((Mdim + mdim + m + k) * sizeof(rsb_coo_idx_t)) / (nnz * el_size)) > 8 * cs)
			sp = RSB_BOOL_TRUE;
		else
			sp = RSB_BOOL_FALSE;
	}
	else
	{
		/* no cache info (FIXME: there should be no section like this one) */
		if (
			Mdim < 8 || mdim < 8 || m < 500 || k < 500 || nnz < 200 * 100)
			sp = RSB_BOOL_FALSE;
		else
			sp = RSB_BOOL_TRUE;
	}
#ifdef RSB_EXPERIMENTAL_ROWS_SUBDIVIDE_TO_CORES_NUM
	/* STILL UNIMPLEMENTED */
#endif /* RSB_EXPERIMENTAL_ROWS_SUBDIVIDE_TO_CORES_NUM */
#if RSB_EXPERIMENTAL_NO_SUBDIVIDE_ON_MIN_NNZ_PER_ROW_OR_COLUMN
	if (1)
	{
		rsb_nnz_idx_t nnzpr;
		if ((flags & RSB_FLAG_WANT_COLUMN_MAJOR_ORDER) != 0)
			nnzpr = nnz / k;
		else
			nnzpr = nnz / m;

		if (nnzpr < RSB_CONST_MIN_NNZ_PER_ROW_OR_COLUMN_PER_SUBMATRIX)
			sp = RSB_BOOL_FALSE;
	}
#endif /* RSB_EXPERIMENTAL_NO_SUBDIVIDE_ON_MIN_NNZ_PER_ROW_OR_COLUMN */
done:
	return sp;
#else /* (RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY == RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE) */
#error "should use a RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE partitioning policy!"
	return RSB_BOOL_FALSE;
#endif /* (RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY == RSB_EXPERIMENTAL_QUAD_DIVISION_POLICY_NAIVE) */
}

static rsb_err_t rsb_do_compute_vertical_split_search_only(
	const rsb_coo_idx_t *RSB_RESTRICT IA, const rsb_coo_idx_t *RSB_RESTRICT JA,
	rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_coo_idx_t m, rsb_coo_idx_t k,
	rsb_coo_idx_t hm, rsb_coo_idx_t hk, rsb_nnz_idx_t nnz,
	const rsb_coo_idx_t *IB, rsb_nnz_idx_t *ulp, rsb_nnz_idx_t *urp, rsb_nnz_idx_t *llp, rsb_nnz_idx_t *lrp)
{
	/**
	\ingroup gr_unfinished


	 */
	rsb_nnz_idx_t ul = 0, ur = 0, ll = 0, lr = 0;
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_nnz_idx_t dnnz = 0 /*,wdnnz = 0,rnnz = 0,hrnnz = 0*/;
	register rsb_coo_idx_t i;
	// rsb_nnz_idx_t nnz0 = 0;
	// rsb_coo_idx_t xroff = 0;

	if (nnz > m || 1)
	// if(nnz>m)
	{
		for (i = roff; RSB_LIKELY(i < roff + m); ++i)
		{
			// offset of line i in the global line pointers array
			rsb_nnz_idx_t nnz0 = IB[i];
			// nnz1..nnz0 are the boundaries of line i
			rsb_nnz_idx_t nnz1 = IB[i + 1];
			rsb_nnz_idx_t nnz2 = 0;
			// check
			assert(nnz0 >= IB[i]);
			assert(nnz1 <= IB[i + 1]);
			// skip line if empty
			if (nnz1 - nnz0 < 1)
				continue;
			// find first element of line i also in the submatrix
			nnz0 += rsb__nnz_split_coo_bsearch(JA + nnz0, coff, nnz1 - nnz0);
			// skip line if empty in the submatrix
			if (nnz1 - nnz0 < 1)
				continue;
			// find the length of the subrow i in the submatrix
			nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, coff + k, nnz1 - nnz0);
			// check
			assert(JA[nnz0 + 0] >= coff);
			// skip line if empty in the submatrix
			if (nnz1 - nnz0 < 1)
				continue;
			nnz2 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, coff + hk, nnz1 - nnz0);
			assert(nnz1 <= IB[i + 1]);
			assert(JA[nnz0 + 0] >= coff);
			assert(JA[nnz1 - 1] < coff + k);
			dnnz += nnz1 - nnz0;
			if (i < roff + hm)
				ul += nnz2 - nnz0,
					ur += nnz1 - nnz2;
			else
				ll += nnz2 - nnz0,
					lr += nnz1 - nnz2;
		}
	}
	else
	{
		// FIXME: UNFINISHED
		rsb_nnz_idx_t nnz0, nnz1, n;
		// RSB_INFO("almost empty matrix !\n");
		for (n = 0; n < nnz; ++n)
		{
			rsb_nnz_idx_t nnz2 = 0;
			i = IA[n];
			nnz0 = IB[i];
			nnz1 = IB[i + 1];
			// ...

			// skip line if empty
			if (nnz1 - nnz0 < 1)
				continue;
			// find first element of line i also in the submatrix
			nnz0 += rsb__nnz_split_coo_bsearch(JA + nnz0, coff, nnz1 - nnz0);
			// skip line if empty in the submatrix
			if (nnz1 - nnz0 < 1)
				continue;
			// find the length of the subrow i in the submatrix
			nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, coff + k, nnz1 - nnz0);
			// check
			assert(JA[nnz0 + 0] >= coff);
			// skip line if empty in the submatrix
			if (nnz1 - nnz0 < 1)
				continue;
			nnz2 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, coff + hk, nnz1 - nnz0);
			assert(nnz1 <= IB[i + 1]);
			assert(JA[nnz0 + 0] >= coff);
			assert(JA[nnz1 - 1] < coff + k);
			dnnz += nnz1 - nnz0;
			if (i < roff + hm)
				ul += nnz2 - nnz0,
					ur += nnz1 - nnz2;
			else
				ll += nnz2 - nnz0,
					lr += nnz1 - nnz2;
		}
	}
	// done:
	*llp = ll;
	*lrp = lr;
	*ulp = ul;
	*urp = ur;
	// err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_fill_early_leaf_matrix(struct rsb_mtx_t *mtxAp, struct rsb_mtx_t *submatrix,
											   const rsb_coo_idx_t *IL, const rsb_coo_idx_t *IR,
											   rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, const rsb_coo_idx_t *VA,
											   // const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_coo_idx_t * VA,
											   rsb_nnz_idx_t snzoff, rsb_nnz_idx_t nnz, rsb_coo_idx_t m, rsb_coo_idx_t k,
											   rsb_coo_idx_t roff, rsb_coo_idx_t coff, rsb_type_t typecode, rsb_flags_t flags)
{
	/**
		\ingroup gr_unfinished
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;

#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
	RSB_INFO("building a very sparse recursive matrix\n");
#endif

	/* no hope for CSR : however, full/half word COO will fit  */
	submatrix->nzoff = snzoff;
	submatrix->bindx = NULL;
	submatrix->bpntr = NULL;
	RSB_DO_FLAG_DEL(submatrix->flags, RSB_FLAG_WANT_BCSS_STORAGE);
	// RSB_ERROR("nnz=%d ,m=%d ! what shall we do ?\n",nnz,m);

	mtxAp->sm[roff ? (coff ? 3 : 2) : (coff ? 1 : 0)] = submatrix;
	submatrix->roff = roff + mtxAp->roff;
	submatrix->coff = coff + mtxAp->coff;
	RSB_DO_ERROR_CUMULATE(errval, rsb__set_init_flags_and_stuff(submatrix, NULL, NULL, m, k, nnz, nnz, nnz, typecode, flags));
	// err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_coo2rec_subdivide_parallel(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t *pinfop, rsb_flags_t flags, rsb_err_t *errvalp, struct rsb_mtx_t **submatricesp, struct rsb_mtx_t *mtxAp, const rsb_nnz_idx_t *IB, const rsb_nnz_idx_t *IX, rsb_coo_idx_t *IT, rsb_coo_idx_t *WA, rsb_submatrix_idx_t cmc, rsb_submatrix_idx_t omc, rsb_submatrix_idx_t tmc, rsb_thread_t wet, rsb_submatrix_idx_t *cmcp)
{
	/*
		TODO: clean this up.
		Note that rsb__set_num_threads is outlawed here.
	 */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t el_size = rsb__sizeof(typecode);
	const rsb_nnz_idx_t ttlnz = nnz;					/* total nnz */
	rsb_nnz_idx_t maxnz = nnz;							/* max encountered nnz for a leaf */
	rsb_submatrix_idx_t stmc = RSB_MIN(tmc, wet);		/* submatrices total count */
	rsb_submatrix_idx_t lmc = 1;						/* leaf matrix count */
	rsb_time_t cpt = RSB_TIME_ZERO, dt = RSB_TIME_ZERO; /* cpt overwrites mtxAp->cpt */
	// const rsb_thread_t tn = rsb_get_num_coo2rec_threads(); /* threads number */
	const rsb_thread_t tn = nnz; // 1 (global) CUDA thread for each element
	const rsb_thread_t mtn = RSB_MAX(1, (rsb_thread_t)(rsb_global_session_handle.subdivision_multiplier * tn));
	rsb_thread_t tnn = 1;									 /* threads number */
	rsb_float_t skew = ((rsb_float_t)(maxnz)) / (nnz / wet); /* if more than one, will limit scaling */
	const long cbs = rsb__get_cache_block_byte_size();
#if defined(RSB_WANT_VERBOSE_SUBDIVISION) && RSB_WANT_VERBOSE_SUBDIVISION != 0
	RSB_INFO("serial substage subdivision of "), RSB_INFO_MATRIX_SUMMARY(mtxAp), RSB_INFO("\n");
#endif
again:
#pragma omp parallel reduction(| \
							   : errval) shared(submatricesp) num_threads(tnn)
{
	rsb_submatrix_idx_t smi = 0;
	struct rsb_mtx_t *submatrix = NULL;
#if RSB_WANT_OMP_RECURSIVE_KERNELS
	rsb_thread_t th_id = omp_get_thread_num();
#else  /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	rsb_thread_t th_id = 0;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#if RSB_WANT_VERBOSE_SUBDIVISION2
	if (th_id == 0)
	{
		RSB_MTXASM_INFO("entering %d threads in subdivision phase\n", tnn);
		RSB_MTXASM_INFO("entering %d threads in subdivision phase\n", tnn);
	}
#endif /* RSB_WANT_VERBOSE_SUBDIVISION2 */
iagain:
#pragma omp critical(rsb_coo2rsbsub_crs)
{
	smi = rsb_do_pick_largest_open_matrix(submatricesp + cmc, omc);
	if (smi != RSB_SUBM_IDX_MARKER)
	{
		smi += cmc;
		submatrix = submatricesp[smi];
		maxnz = submatrix->nnz;
		/* RSB_ASSERT(nnz>=wet); */
		skew = ((rsb_float_t)(maxnz)) / ((rsb_float_t)(nnz / wet));
		RSB_ASSERT(!isinf(skew));
		omc--;
		if (smi != cmc)
		{
			RSB_SWAP(struct rsb_mtx_t *, submatricesp[smi], submatricesp[cmc]);
		}
		++cmc;
	}
	else
	{
		submatrix = NULL;
	}

#if defined(RSB_WANT_VERBOSE_SUBDIVISION) && RSB_WANT_VERBOSE_SUBDIVISION != 0
	if (submatrix)
		RSB_INFO("subdividing "), RSB_INFO_MATRIX_SUMMARY(submatrix), RSB_INFO(" (open:%d,closed:%d) for thread %d\n", omc, cmc, th_id);
	else
		RSB_INFO("no available submatrix (open:%d,closed:%d) for thread %d/%d\n", omc, cmc, th_id, tnn);
#endif
}
	if ((smi) != RSB_SUBM_IDX_MARKER)
	{
		{
			const rsb_coo_idx_t k = submatrix->nc;
			const rsb_coo_idx_t m = submatrix->nr;
			const rsb_coo_idx_t hk = RSB_MIDDLE(k);
			const rsb_coo_idx_t hm = RSB_MIDDLE(m);
			rsb_nnz_idx_t ul = 0, ur = 0, ll = 0, lr = 0;
			const rsb_nnz_idx_t nnz = submatrix->nnz;
			const rsb_coo_idx_t roff = submatrix->roff;
			const rsb_coo_idx_t coff = submatrix->coff;
			const rsb_flags_t flags = submatrix->flags;
			const rsb_nnz_idx_t nzoff = submatrix->nzoff;
			rsb_bool_t sqs = RSB_BOOL_FALSE; /* should quad subdivide */
			rsb_submatrix_idx_t smc = 0;	 /* submatrices count */

#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
			RSB_INFO("cmc:%d omc:%d smi:%d tmc=%d stmc=%d th_id=%d\n", cmc, omc, smi, tmc, stmc, th_id);
#endif

			/* too few nonzeros for recursion (TODO: may change in the future) */
			if (RSB_DO_TOOFEWNNZFORRCSR(nnz, m))
#if RSB_WANT_SUBDIVISION_FIXES_20101120
				if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_WANT_COO_STORAGE))
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */
				{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
					RSB_INFO("matrix too sparse for RCSR: rejoining\n");
#endif
					sqs = RSB_BOOL_FALSE;
					goto nosqstest;
				}

			/* decide if the matrix is worth subdividing further (soft) */
			sqs = rsb__should_recursively_partition_matrix(0, 0, m, k, 0, 0, nnz, m, k, roff, coff, flags, el_size, mtn);
#if RSB_WANT_SUBDIVISION_FIXES_20101120
			if (nnz < RSB_RECURSION_MIN_NNZ || m < RSB_RECURSION_MIN_DIM || k < RSB_RECURSION_MIN_DIM || !RSB_DO_FLAG_HAS(flags, RSB_FLAG_QUAD_PARTITIONING))
			{
				sqs = RSB_BOOL_FALSE; /* a hard condition */
				goto nosqstest;
			}
			else if (cmc + omc < tmc)
				if (skew > RSB_SUBDIVISION_SKEW_MAX)
					sqs = RSB_BOOL_TRUE; /* a soft condition */
#endif									 /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

			if (!sqs)
				if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS))
					if (wet > lmc)
						sqs = RSB_BOOL_TRUE;

			if (sqs)
			{
				rsb_bool_t awfcsr = RSB_BOOL_FALSE; /* all of the matrices will fit csr ? */
#if RSB_WANT_SUBDIVISION_FIXES_20101120
				rsb_nnz_idx_t mqnnz = RSB_MAX(RSB_MAX(ul, ur), RSB_MAX(lr, ll));
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

				/* compute the split vector */
				dt = -rsb_time();
				if ((errval = rsb_do_compute_vertical_split_search_only(IA, JA, roff, coff, m, k, hm, hk, nnz, IB, &ul, &ur, &ll, &lr)) != RSB_ERR_NO_ERROR)
					; /* goto err; */
				dt += rsb_time();
				cpt += dt;
				// assert(IR);
				awfcsr = ((ul > 0 && RSB_DO_TOOFEWNNZFORCSR(ul, hm)) || (ur > 0 && RSB_DO_TOOFEWNNZFORCSR(ur, hm)) || (lr > 0 && RSB_DO_TOOFEWNNZFORCSR(lr, m - hm)) || (ll > 0 && RSB_DO_TOOFEWNNZFORCSR(ll, m - hm))) ? RSB_BOOL_TRUE : RSB_BOOL_FALSE;

				if (awfcsr) /* FIXME: misleading naming ! */
				{
					/* if some leaf won't fit in CSR, we don't split anymore */
					if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_WANT_COO_STORAGE))
						sqs = RSB_BOOL_TRUE;
					else
						sqs = RSB_BOOL_FALSE;
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
					RSB_INFO("no space for conversion of some leaf: rejoining ? %d\n", !sqs);
#endif
				}

#if RSB_WANT_SUBDIVISION_FIXES_20101120
				/* an alternative would be to place this test in the branch above*/
				if (mqnnz > RSB_MAX_QUADRANTS_UNBALANCE * (nnz - mqnnz) &&
					el_size * nnz < cbs &&
					nnz < (ttlnz / wet))
					sqs = RSB_BOOL_FALSE;
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

				/* how many submatrices out of four ? */
				smc = (ul ? 1 : 0) + (ur ? 1 : 0) + (ll ? 1 : 0) + (lr ? 1 : 0);
				if (cmc + omc + smc > tmc)
				{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
					RSB_INFO("too many submatrices (%d+%d>%d: rejoining\n", cmc + omc, smc, tmc);
#endif
					sqs = RSB_BOOL_FALSE;
					goto nosqstest;
				}

#if !RSB_WANT_SUBDIVISION_FIXES_20101120
#ifdef RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES
				if (RSB_DO_FLAG_HAS(flags, RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES))
					if (wet < lmc - 1)
					{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
						RSB_INFO("RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES: rejoining\n");
#endif
						sqs = RSB_BOOL_FALSE;
					}
#endif /* RSB_FLAG_EXPERIMENTAL_NO_MICRO_LEAVES */
#endif /* RSB_WANT_SUBDIVISION_FIXES_20101120 */

#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_ERROR("splitting %d/%d -> %d/%d %d/%d %d/%d %d/%d sqs? %d\n", nnz, m, ul, hm, ur, hm, ll, m - hm, lr, m - hm, sqs);
#endif
				if (ul + ur + ll + lr != nnz)
				{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
					RSB_ERROR("%d ?= %d + %d + %d + %d = %d\n", nnz, ul, ur, ll, lr, ul + ur + ll + lr);
#endif
					RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
				}
			}
		nosqstest:
			if (sqs)
			{
				/* should quad-subdivide. let's take care of indices. */
				rsb_nnz_idx_t snzoff = nzoff;
				rsb_submatrix_idx_t smci = 0;
				rsb_submatrix_idx_t smco = 0;
				struct rsb_mtx_t *isms[4] = {NULL, NULL, NULL, NULL};

				/*
				the index arrays are copied/linked into the quadrants
				some quadrants may seem ready for recursion, but they not result as such later on.
				they will be made leaf later on, if necessary.
				...
				*/
				assert(ur >= 0 && ul >= 0 && lr >= 0 && ll >= 0);

#pragma omp critical(rsb_coo2rsbsub_crs)
				{
					if (cmc + omc + smc + RSB_SUBDIVISION_BUG_EXTRA > tmc)
					{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
						RSB_INFO("too many submatrices (%d+%d>%d): rejoining\n", cmc + omc, smc, tmc);
#endif
						sqs = RSB_BOOL_FALSE;
					}
					else
					{
						lmc += smc;
						lmc -= 1;
						smco = cmc + omc;
						snzoff = nzoff;
						if (ul)
						{
							isms[0] = submatricesp[smco + smci];
							submatricesp[smco + smci] = NULL;
							smci++;
						}
						if (ur)
						{
							isms[1] = submatricesp[smco + smci];
							submatricesp[smco + smci] = NULL;
							smci++;
						}
						if (ll)
						{
							isms[2] = submatricesp[smco + smci];
							submatricesp[smco + smci] = NULL;
							smci++;
						}
						if (lr)
						{
							isms[3] = submatricesp[smco + smci];
							submatricesp[smco + smci] = NULL;
							smci++;
						}
						smci = 0;
						omc += smc;
					}
					if (sqs)
					{
						if (ul)
							RSB_DO_ERROR_CUMULATE(errval, rsb_do_fill_early_leaf_matrix(submatrix, isms[0], NULL, NULL, IA, JA, (const rsb_coo_idx_t *)VA, snzoff, ul, hm, hk, 0, 0, typecode, flags)), snzoff += ul, ++smci;
						if (ur)
							RSB_DO_ERROR_CUMULATE(errval, rsb_do_fill_early_leaf_matrix(submatrix, isms[1], NULL, NULL, IA, JA, (const rsb_coo_idx_t *)VA, snzoff, ur, hm, k - hk, 0, hk, typecode, flags)), snzoff += ur, ++smci;
						if (ll)
							RSB_DO_ERROR_CUMULATE(errval, rsb_do_fill_early_leaf_matrix(submatrix, isms[2], NULL, NULL, IA, JA, (const rsb_coo_idx_t *)VA, snzoff, ll, m - hm, hk, hm, 0, typecode, flags)), snzoff += ll, ++smci;
						if (lr)
							RSB_DO_ERROR_CUMULATE(errval, rsb_do_fill_early_leaf_matrix(submatrix, isms[3], NULL, NULL, IA, JA, (const rsb_coo_idx_t *)VA, snzoff, lr, m - hm, k - hk, hm, hk, typecode, flags)), snzoff += lr, ++smci;
						RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
					}

					if (sqs)
					{
						smci = 0;
						if (ul)
						{
							submatricesp[smco + smci] = isms[0];
							smci++;
						}
						if (ur)
						{
							submatricesp[smco + smci] = isms[1];
							smci++;
						}
						if (ll)
						{
							submatricesp[smco + smci] = isms[2];
							smci++;
						}
						if (lr)
						{
							submatricesp[smco + smci] = isms[3];
							smci++;
						}
					}

					if (sqs)
					{
						if (snzoff - nzoff != nnz)
						{
							/* is this a partition ? */
							RSB_ERROR("%d - %d != %d ?= %d + %d + %d + %d = %d\n", snzoff, nzoff, nnz, ul, ur, ll, lr, ul + ur + ll + lr);
							RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
						}
						if (RSB_SOME_ERROR(errval))
						{
							RSB_ERROR(RSB_ERRM_ES); /* goto err; */
						}
						RSB_DO_FLAG_ADD(submatrix->flags, RSB_FLAG_QUAD_PARTITIONING);
						RSB_DO_FLAG_DEL(submatrix->flags, RSB_FLAG_NON_ROOT_MATRIX);
						submatrix->bindx = NULL;
						submatrix->bpntr = NULL;
						submatrix->indptr = NULL;
					}
				}
			}
			if (!sqs)
			{
				RSB_DO_FLAG_SUBST(submatrix->flags, RSB_FLAG_QUAD_PARTITIONING, RSB_FLAG_NON_ROOT_MATRIX);
				/* selecting a format and declaring as leaf */
				if (!RSB_DO_TOOFEWNNZFORCSR(nnz, m) /*&& IR && IL*/)
				{
					/*				RSB_INFO("CSR -> COO ?\n"); */
					if (RSB_DO_FLAG_HAS(submatrix->flags, RSB_FLAG_WANT_BCSS_STORAGE))
						RSB_DO_FLAG_DEL(submatrix->flags, RSB_FLAG_WANT_COO_STORAGE);
					if ((errval = rsb__do_set_init_storage_flags(submatrix, submatrix->flags)) != RSB_ERR_NO_ERROR)
						; /* goto err; */
				}
				else
				{
					/*				RSB_INFO("COO !\n"); */
					rsb_flags_t sflags = flags;
					RSB_DO_FLAG_SUBST(sflags, RSB_FLAG_WANT_BCSS_STORAGE, RSB_FLAG_WANT_COO_STORAGE);
					if ((errval = rsb__do_set_init_storage_flags(submatrix, sflags)) != RSB_ERR_NO_ERROR)
						; /* goto err; */
				}
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_INFO("freezing %d ", smi + 1),
					RSB_INFO_MATRIX_SUMMARY(submatrix),
					RSB_INFO("\n");
#endif
			}
			/* matrix is declared as 'closed'.
			   sorting, in a way smi will point to the biggest open mtxAp, which will be picked up next */
#if 0
		qsort(submatricesp+cmc,(size_t)(omc),sizeof(struct rsb_mtx_t*),& rsb_compar_rcsr_matrix_regarding_nnz);
		RSB_SWAP(struct rsb_mtx_t *,submatricesp[smi],submatricesp[cmc-1]);
#endif
		}
		/*smi = cmc+omc-1; */
	}
	if (omc > 0 && cmc + omc < stmc)
		goto iagain;
#if RSB_WANT_VERBOSE_SUBDIVISION2
	if (th_id == 0)
	{
		RSB_MTXASM_INFO("thread %d:terminating subdivision", th_id);
		if (omc == 0)
		{
			RSB_MTXASM_INFO(", no more open matrices");
		}
		RSB_MTXASM_INFO("(closed %d= %d nodes + %d leaves, out of %d available)", cmc, cmc - lmc, lmc, tmc);
		RSB_MTXASM_INFO(",(maxnz=%d,skew=%g)", maxnz, skew);
		if (cmc + omc >= stmc)
		{
			RSB_MTXASM_INFO(", no room left for submatrices");
		}
		RSB_MTXASM_INFO(".\n");
	}
#endif
} /* parallel */

	if (RSB_SOME_ERROR(errval))
		goto err;

#pragma omp barrier
	if (stmc != tmc)
	{
		stmc = tmc;
#if defined(RSB_WANT_VERBOSE_SUBDIVISION) && RSB_WANT_VERBOSE_SUBDIVISION != 0
		RSB_INFO("parallel substage subdivision of "), RSB_INFO_MATRIX_SUMMARY(mtxAp), RSB_INFO("\n");
#endif
		tnn = tn;
		goto again;
	}
	else
	{

#if defined(RSB_WANT_VERBOSE_SUBDIVISION) && RSB_WANT_VERBOSE_SUBDIVISION != 0
		RSB_INFO("parallel substage subdivision of "), RSB_INFO_MATRIX_SUMMARY(mtxAp), RSB_INFO(" not required\n");
#endif
	}
	{
#if defined(RSB_WANT_VERBOSE_SUBDIVISION) && RSB_WANT_VERBOSE_SUBDIVISION != 0
		RSB_INFO("subdivision of "), RSB_INFO_MATRIX_SUMMARY(mtxAp), RSB_INFO("complete \n");
#endif
	}
	mtxAp->cpt = cpt;

	*cmcp = cmc;
err:
	RSB_DO_ERR_RETURN(errval)
}

static int rsb__compar_rcsr_matrix_leftmost_first(const void *ap, const void *bp)
{
	struct rsb_mtx_t *a = *(struct rsb_mtx_t **)ap;
	struct rsb_mtx_t *b = *(struct rsb_mtx_t **)bp;
	const rsb_bool_t at = !RSB_DO_FLAG_HAS(a->flags, RSB_FLAG_QUAD_PARTITIONING);
	const rsb_bool_t bt = !RSB_DO_FLAG_HAS(b->flags, RSB_FLAG_QUAD_PARTITIONING);
	int ss = 1; /* should swap results ? */

	if (at && !bt)
		return 1;
	if (!at && bt)
		return -1;

	if (a->coff < b->coff)
	{
		RSB_SWAP(struct rsb_mtx_t *, a, b);
		ss = -1; /* should swap results ! */
	}

	return (a->coff == b->coff) ? (a->roff > b->roff ? 1 : (a->roff < b->roff ? -1 : 0)) : 1 * ss;
}

static rsb_err_t rsb_do_coo2rec_shuffle(void *VA, rsb_coo_idx_t *IA, rsb_coo_idx_t *JA, rsb_coo_idx_t m, rsb_coo_idx_t k, rsb_nnz_idx_t nnz, rsb_type_t typecode, const struct rsb_mtx_partitioning_info_t *pinfop, rsb_flags_t flags, rsb_err_t *errvalp, struct rsb_mtx_t **submatricesp, struct rsb_mtx_t *mtxAp, const rsb_nnz_idx_t *IB, rsb_coo_idx_t *WA, rsb_submatrix_idx_t cmc)
{
	rsb_nnz_idx_t tdnnz = 0;
	rsb_submatrix_idx_t smi = 0; /* max matrix count, done matrix count, submatrix index */
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	size_t el_size = rsb__sizeof(typecode);
#if RSB_WANT_VERBOSE_TIMINGS
	rsb_time_t pmt = RSB_TIME_ZERO;
#endif /* RSB_WANT_VERBOSE_TIMINGS */

	if (!VA || cmc == 1)
	{
		mtxAp->VA = VA;
		goto no_va_cp;
	}

// the following is a highly parallel phase
#pragma omp parallel for schedule(static, 1) reduction(| \
													   : errval) shared(tdnnz, submatricesp, IB) RSB_NTC
	for (smi = 0; smi < cmc; ++smi)
	{
		struct rsb_mtx_t *submatrix = submatricesp[smi];
#if RSB_OBSOLETE_QUARANTINE
		const rsb_coo_idx_t *IL = submatrix->bindx;
		const rsb_coo_idx_t *IR = submatrix->bpntr;
#endif /* RSB_OBSOLETE_QUARANTINE */
		rsb_nnz_idx_t dnnz = 0;
		rsb_coo_idx_t i;
		//		if(RSB_C2R_IF_VERBOSE)
		//		RSB_INFO("%d -> %d\n",smi,omp_get_thread_num());
		if (!RSB_DO_FLAG_HAS((submatrix->flags), RSB_FLAG_QUAD_PARTITIONING))
		{
#if RSB_OBSOLETE_QUARANTINE
			if (RSB_DO_ENOUGHNNZFORINDEXBASEDBUILD(submatrix) && !rsb__is_coo_matrix(submatrix) && IL && IR)
			{
				assert(IL);
				assert(IR);
				dnnz = rsb_do_copy_submatrix_coa(submatrix, VA, WA, IL, IR, el_size, 0, 0, submatrix->nr);
#if 0
				for(i=submatrix->roff;RSB_LIKELY(i<submatrix->roff+submatrix->nr);++i)
				{
					rsb_nnz_idx_t nnz1 = IR[i-submatrix->roff];
					rsb_nnz_idx_t nnz0 = IL[i-submatrix->roff];
					assert(IL[i-submatrix->roff]>=IB[i]);
					assert(IL[i-submatrix->roff]<=IB[i+1]);
					assert(IL[i-submatrix->roff+1]>=IB[i+1]);
					assert(nnz0>=IL[i-submatrix->roff]);
					assert(nnz1<=IL[i-submatrix->roff+1]);
				       	assert(nnz1<=IB[i+1]);
					if(IB[i]==IB[i+1])continue;
					assert(nnz1>=nnz0);
					if(nnz1==nnz0)continue;
					assert(JA[nnz0+0]>=submatrix->coff);
					assert(JA[nnz1-1]< submatrix->coff+submatrix->nc);
				}
#endif
			}
			else
#endif /* RSB_OBSOLETE_QUARANTINE */
				if (!RSB_DO_TOOFEWNNZFORCSR(submatrix->nnz, submatrix->nr))
				{
					// rsb_coo_idx_t*IL = submatrix->bindx;
					for (i = submatrix->roff; RSB_LIKELY(i < submatrix->roff + submatrix->nr); ++i)
					{
						// rsb_nnz_idx_t fel;
						//  offset of line i in the global line pointers array
						rsb_nnz_idx_t nnz0 = IB[i];
						// nnz1..nnz0 are the boundaries of line i
						rsb_nnz_idx_t nnz1 = IB[i + 1];
						// check
						assert(nnz0 >= IB[i]);
						assert(nnz1 <= IB[i + 1]);
						// skip line if empty
						if (nnz1 - nnz0 < 1)
							continue;
						// find first element of line i also in the submatrix
						nnz0 += rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff, nnz1 - nnz0);
						// skip line if empty in the submatrix
						if (nnz1 - nnz0 < 1)
							continue;
						// find the length of the subrow i in the submatrix
						nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff + submatrix->nc, nnz1 - nnz0);
						// check
						assert(JA[nnz0 + 0] >= submatrix->coff);
						// skip line if empty in the submatrix
						if (nnz1 - nnz0 < 1)
							continue;
						// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
						//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
						// checks
						//					assert(IL[i-submatrix->roff]>=IB[i]);
						//					assert(IL[i-submatrix->roff]<=IB[i+1]);
						//					assert(IL[i-submatrix->roff+1]>=IB[i+1]);
						//					assert(nnz0>=IL[i-submatrix->roff]);
						//					assert(nnz1<=IL[i-submatrix->roff+1]);
						assert(nnz1 <= IB[i + 1]);
						assert(JA[nnz0 + 0] >= submatrix->coff);
						assert(JA[nnz1 - 1] < submatrix->coff + submatrix->nc);
						// perform the copy
						RSB_A_MEMCPY_SMALL(WA, VA, submatrix->nzoff + dnnz, nnz0, nnz1 - nnz0, el_size);
						// RSB_COA_MEMCPY(WA,JA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0);
						//  update the actual offset in the destination array
						dnnz += nnz1 - nnz0;
					}
				}
				else
				{
					// rsb_coo_idx_t*IL = submatrix->bindx;
					for (i = submatrix->roff; RSB_LIKELY(i < submatrix->roff + submatrix->nr); ++i)
					{
						// rsb_nnz_idx_t fel;
						//  offset of line i in the global line pointers array
						rsb_nnz_idx_t nnz0 = IB[i];
						// nnz1..nnz0 are the boundaries of line i
						rsb_nnz_idx_t nnz1 = IB[i + 1];
						// check
						assert(nnz0 >= IB[i]);
						assert(nnz1 <= IB[i + 1]);
						// skip line if empty
						if (nnz1 - nnz0 < 1)
							continue;
						// find first element of line i also in the submatrix
						nnz0 += rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff, nnz1 - nnz0);
						// skip line if empty in the submatrix
						if (nnz1 - nnz0 < 1)
							continue;
						// find the length of the subrow i in the submatrix
						nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff + submatrix->nc, nnz1 - nnz0);
						// check
						assert(JA[nnz0 + 0] >= submatrix->coff);
						// skip line if empty in the submatrix
						if (nnz1 - nnz0 < 1)
							continue;
						// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
						//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
						// checks
						//					assert(IL[i-submatrix->roff]>=IB[i]);
						//					assert(IL[i-submatrix->roff]<=IB[i+1]);
						//					assert(IL[i-submatrix->roff+1]>=IB[i+1]);
						//					assert(nnz0>=IL[i-submatrix->roff]);
						//					assert(nnz1<=IL[i-submatrix->roff+1]);
						assert(nnz1 <= IB[i + 1]);
						assert(JA[nnz0 + 0] >= submatrix->coff);
						assert(JA[nnz1 - 1] < submatrix->coff + submatrix->nc);
						// perform the copy
						RSB_A_MEMCPY_SMALL(WA, VA, submatrix->nzoff + dnnz, nnz0, nnz1 - nnz0, el_size);
						// RSB_COA_MEMCPY(WA,JA,submatrix->nzoff+dnnz,nnz0,nnz1-nnz0);
						//  update the actual offset in the destination array
						dnnz += nnz1 - nnz0;
					}
				}
			if (dnnz != submatrix->nnz)
			{
				RSB_ERROR("@%d,%d: found %d, should have found %d\n",
						  submatrix->roff, submatrix->coff, dnnz, submatrix->nnz);
				RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
			}
#pragma omp critical(rsb_coo2rsb_nzinc_crs)
			{
				tdnnz += dnnz;
			}
		}
	}
	// gerr:
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS
#if !defined(__xlC__)
/* FIXME: xlc does not allow this, but we have experienced problems, without */
#pragma omp barrier
#endif /* __xlC__ */
	if (RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err, RSB_ERRM_NL);
	}

	if (tdnnz != nnz)
	{
		RSB_ERROR("found %d, should have found %d\n", tdnnz, nnz);
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}

#if RSB_WANT_VERBOSE_TIMINGS
	pmt -= rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */
	RSB_A_MEMCPY_parallel(VA, WA, 0, 0, nnz, el_size);
#if RSB_WANT_VERBOSE_TIMINGS
	pmt += rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */
no_va_cp:

	tdnnz = 0;
#pragma omp parallel for schedule(static, 1) reduction(| \
													   : errval) shared(tdnnz, submatricesp, IB) RSB_NTC
	for (smi = 0; smi < cmc; ++smi)
	{
		struct rsb_mtx_t *submatrix = submatricesp[smi];
		if (!RSB_DO_FLAG_HAS(submatrix->flags, RSB_FLAG_QUAD_PARTITIONING))
		{
			rsb_coo_idx_t oll, nll;
			rsb_coo_idx_t i;
			const rsb_coo_idx_t *IR = submatrix->bpntr;
			const rsb_coo_idx_t *IL = submatrix->bindx;
			rsb_nnz_idx_t dnnz = 0;

			// if(!RSB_DO_TOOFEWNNZFORRCSR(submatrix->nnz,submatrix->nr))
			if (RSB_DO_ENOUGHNNZFORINDEXBASEDBUILD(submatrix) && !rsb__is_coo_matrix(submatrix) && IL && IR)
			{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_INFO("CSR:%zd/%zd:%zd..%zd\n", (rsb_printf_int_t)smi, (rsb_printf_int_t)cmc, (rsb_printf_int_t)submatrix->nzoff, (rsb_printf_int_t)(submatrix->nzoff + submatrix->nnz));
#endif
				submatrix->bpntr = IA + submatrix->nzoff;
				assert(IL);
				assert(IR);
				assert(IR < IA + submatrix->nzoff + submatrix->nnz);
				assert(IL < IA + submatrix->nzoff + submatrix->nnz);
				assert(IR >= IA + submatrix->nzoff);
				assert(IL >= IA + submatrix->nzoff);
				for (dnnz = 0, i = 0; RSB_LIKELY(i < submatrix->nr); dnnz += IR[i] - IL[i], ++i)
					RSB_COA_MEMCPY_SMALL(WA, JA, submatrix->nzoff + dnnz, IL[i], IR[i] - IL[i]);

				//				RSB_INFO("%d x %d (%d) @ %d, %d (rcsr)\n",submatrix->nr,submatrix->nc,submatrix->nnz,submatrix->roff,submatrix->coff);
				if (dnnz != submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d, should have found %d\n",
							  submatrix->roff, submatrix->coff, dnnz, submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
				}
				assert(IL);
				assert(IR);

				oll = IR[0] - IL[0];
				submatrix->bpntr[0] = 0;
				for (i = 1; RSB_LIKELY(i < submatrix->nr); ++i)
				{
					nll = IR[i] - IL[i];
					submatrix->bpntr[i] = submatrix->bpntr[i - 1] + oll;
					oll = nll;
				}
				submatrix->bpntr[submatrix->nr] = submatrix->bpntr[submatrix->nr - 1] + oll;
				if (submatrix->bpntr[submatrix->nr] != submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d, should have found %d\n",
							  submatrix->roff, submatrix->coff, submatrix->bpntr[submatrix->nr], submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
				}
			}
			else if (!RSB_DO_TOOFEWNNZFORCSR(submatrix->nnz, submatrix->nr) && !rsb__is_coo_matrix(submatrix))
			{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_INFO("CSR:%zd/%zd:%zd..%zd\n", (rsb_printf_int_t)smi, (rsb_printf_int_t)cmc, (rsb_printf_int_t)submatrix->nzoff, (rsb_printf_int_t)(submatrix->nzoff + submatrix->nnz));
#endif
				oll = 0;
				submatrix->bpntr = IA + submatrix->nzoff;
				submatrix->bpntr[0] = 0;
				//				RSB_INFO("%d x %d (%d) @ %d, %d (rcsr)\n",submatrix->nr,submatrix->nc,submatrix->nnz,submatrix->roff,submatrix->coff);
				for (i = submatrix->roff; RSB_LIKELY(i < submatrix->roff + submatrix->nr); ++i)
				{
					// rsb_nnz_idx_t fel;
					//  offset of line i in the global line pointers array
					rsb_nnz_idx_t nnz0 = IB[i];
					// nnz1..nnz0 are the boundaries of line i
					rsb_nnz_idx_t nnz1 = IB[i + 1];
					// check
					assert(nnz0 >= IB[i]);
					assert(nnz1 <= IB[i + 1]);
					// skip line if empty
					if (nnz1 - nnz0 < 1)
						goto is_empty_subrow;
					// find first element of line i also in the submatrix
					nnz0 += rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff, nnz1 - nnz0);
					// skip line if empty in the submatrix
					if (nnz1 - nnz0 < 1)
						goto is_empty_subrow;
					// find the length of the subrow i in the submatrix
					nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff + submatrix->nc, nnz1 - nnz0);
					// check
					assert(JA[nnz0 + 0] >= submatrix->coff);
					// skip line if empty in the submatrix
					if (nnz1 - nnz0 < 1)
						goto is_empty_subrow;
					// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
					//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
					// checks
					assert(nnz1 <= IB[i + 1]);
					assert(JA[nnz0 + 0] >= submatrix->coff);
					assert(JA[nnz1 - 1] < submatrix->coff + submatrix->nc);
					// convert row indices
					nll = nnz1 - nnz0;
					if (i > submatrix->roff)
						submatrix->bpntr[i - submatrix->roff] = submatrix->bpntr[i - submatrix->roff - 1] + oll;
					oll = nll;
					// perform the copy
					RSB_COA_MEMCPY_SMALL(WA, JA, submatrix->nzoff + dnnz, nnz0, nnz1 - nnz0);
					// update the actual offset in the destination array
					dnnz += nnz1 - nnz0;
					continue;
				is_empty_subrow:
					// convert row indices
					nll = 0;
					if (RSB_LIKELY(i > submatrix->roff))
						submatrix->bpntr[i - submatrix->roff] = submatrix->bpntr[i - submatrix->roff - 1] + oll;
					oll = nll;
				}
				submatrix->bpntr[submatrix->nr] = submatrix->bpntr[submatrix->nr - 1] + oll;
				if (dnnz != submatrix->nnz || submatrix->bpntr[submatrix->nr] != submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d, and %d; should have found %d\n",
							  submatrix->roff, submatrix->coff,
							  dnnz, submatrix->bpntr[submatrix->nr], submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
				}
			}
			else
			{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_INFO("COO:%zd/%zd:%zd..%zd\n", (rsb_printf_int_t)smi, (rsb_printf_int_t)cmc, (rsb_printf_int_t)submatrix->nzoff, (rsb_printf_int_t)(submatrix->nzoff + submatrix->nnz));
#endif
				oll = 0;
				submatrix->bpntr = IA + submatrix->nzoff;
				submatrix->bpntr[0] = 0;
				for (i = submatrix->roff; RSB_LIKELY(i < submatrix->roff + submatrix->nr); ++i)
				{
					// rsb_nnz_idx_t fel;
					//  offset of line i in the global line pointers array
					rsb_nnz_idx_t nnz0 = IB[i];
					// nnz1..nnz0 are the boundaries of line i
					rsb_nnz_idx_t nnz1 = IB[i + 1];
					// check
					assert(nnz0 >= IB[i]);
					assert(nnz1 <= IB[i + 1]);
					// skip line if empty
					if (nnz1 - nnz0 < 1)
						continue;
					// find first element of line i also in the submatrix
					nnz0 += rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff, nnz1 - nnz0);
					// skip line if empty in the submatrix
					if (nnz1 - nnz0 < 1)
						continue;
					// find the length of the subrow i in the submatrix
					nnz1 = nnz0 + rsb__nnz_split_coo_bsearch(JA + nnz0, submatrix->coff + submatrix->nc, nnz1 - nnz0);
					// check
					assert(JA[nnz0 + 0] >= submatrix->coff);
					// skip line if empty in the submatrix
					if (nnz1 - nnz0 < 1)
						continue;
					// nnz1 .. nnz0 contain nonempty subrow i in the submatrix
					//					RSB_INFO("i:%d, %d..%d -> %d\n",i,nnz0,nnz1-1,submatrix->nzoff+dnnz);
					// checks
					assert(nnz1 <= IB[i + 1]);
					assert(JA[nnz0 + 0] >= submatrix->coff);
					assert(JA[nnz1 - 1] < submatrix->coff + submatrix->nc);
					// convert row indices
					// perform the copy
					RSB_COA_MEMCPY_SMALL(WA, JA, submatrix->nzoff + dnnz, nnz0, nnz1 - nnz0);
					rsb__util_coo_array_set(IA + submatrix->nzoff + dnnz, nnz1 - nnz0, i - submatrix->roff);
					// update the actual offset in the destination array
					dnnz += nnz1 - nnz0;
				}
				if (dnnz != submatrix->nnz)
				{
					RSB_ERROR("@%d,%d: found %d; should have found %d\n",
							  submatrix->roff, submatrix->coff, dnnz, submatrix->nnz);
					RSB_DO_ERROR_CUMULATE(errval, RSB_ERR_INTERNAL_ERROR);
				}
			}
			rsb__util_coo_array_sub(WA + submatrix->nzoff, dnnz, submatrix->coff);
#pragma omp critical(rsb_coo2rsb_nzinc_crs)
			{
				tdnnz += dnnz;
			}
		}
	}
#pragma omp barrier
#if RSB_WANT_VERBOSE_TIMINGS
	pmt -= rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */

	RSB_COA_MEMCPY_parallel(JA, WA, 0, 0, nnz);

#if RSB_WANT_VERBOSE_TIMINGS
	pmt += rsb_time();
#endif /* RSB_WANT_VERBOSE_TIMINGS */

	if (RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err, RSB_ERRM_NL);
	}

	if (tdnnz != nnz)
	{
		RSB_ERROR("found %d, should have found %d\n", tdnnz, nnz);
		errval = RSB_ERR_INTERNAL_ERROR;
		goto err;
	}
err:
	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_switch_fresh_terminal_matrix_to_halfword_storages(struct rsb_mtx_t *mtxAp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if (rsb__is_recursive_matrix(mtxAp->flags))
	{
		RSB_ERROR(RSB_ERRM_ES);
		return RSB_ERR_BADARGS;
	}
	else
	{
		if (RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_WANT_COO_STORAGE))
		{
			if (rsb__do_is_candidate_for_halfword_coo(mtxAp))
			{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_INFO("to halfword COO:"), RSB_INFO_MATRIX_SUMMARY(mtxAp), RSB_INFO("\n");
#endif
				RSB_DO_ERROR_CUMULATE(errval, rsb__do_switch_to_halfword_coo(mtxAp));
			}
			else
				RSB_DO_FLAG_DEL(mtxAp->flags, RSB_FLAG_USE_HALFWORD_INDICES);
		}
		else

#if RSB_WANT_FIX_BUG_DISCOVERED_20121210
			if (!RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_WANT_COO_STORAGE))
#else  /* RSB_WANT_FIX_BUG_DISCOVERED_20121210 */
			if (RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_WANT_COO_STORAGE))
#endif /* RSB_WANT_FIX_BUG_DISCOVERED_20121210 */
		{
			if (RSB_SOME_ERROR(rsb__do_is_candidate_for_halfword_csr(mtxAp)))
			{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
				RSB_INFO("to halfword CSR:"), RSB_INFO_MATRIX_SUMMARY(mtxAp), RSB_INFO("\n");
#endif
				RSB_DO_ERROR_CUMULATE(errval, rsb__do_switch_to_halfword_csr(mtxAp));
			}
			else
				RSB_DO_FLAG_DEL(mtxAp->flags, RSB_FLAG_USE_HALFWORD_INDICES);
		}
		else
			// if(!rsb__is_root_matrix(mtxAp) || rsb__is_terminal_recursive_matrix(mtxAp)) /* root recursive or root nonrec. */
			//	RSB_DO_FLAG_DEL(mtxAp->flags,RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES);
			// else
			RSB_DO_FLAG_DEL(mtxAp->flags, RSB_FLAG_USE_HALFWORD_INDICES);
		; /* for root matrices, we keep the flags, because some of the leaves MAY have it */
	}

	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(struct rsb_mtx_t *mtxAp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if (!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		return RSB_ERR_BADARGS;
	}

	if (rsb__is_recursive_matrix(mtxAp->flags))
	{
		rsb_submatrix_idx_t i, j;
		struct rsb_mtx_t *submatrix;
		RSB_SUBMATRIX_FOREACH(mtxAp, submatrix, i, j)
		if (submatrix)
			RSB_DO_ERROR_CUMULATE(errval, rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(submatrix));
	}
	else
		errval = rsb_do_switch_fresh_terminal_matrix_to_halfword_storages(mtxAp);

	RSB_DO_ERR_RETURN(errval)
}

static rsb_err_t rsb_do_compute_bounded_boxes(struct rsb_mtx_t *mtxAp)
{
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_submatrix_idx_t smi = 0;
	rsb_bool_t want_really = 0;
#if RSB_WANT_BOUNDED_BOXES
	want_really = (rsb_global_session_handle.want_bounded_box != 0);
#else  /* RSB_WANT_BOUNDED_BOXES */
	mtxAp->broff = roff;
	mtxAp->bcoff = coff;
	mtxAp->bm = m;
	mtxAp->bk = k;
	goto err;
#endif /* RSB_WANT_BOUNDED_BOXES */

	if (!mtxAp)
	{
		RSB_ERROR(RSB_ERRM_E_MTXAP);
		return RSB_ERR_BADARGS;
	}

	if (mtxAp->nnz == 0)
		goto err;

	if (want_really)
	{
		if (rsb__is_terminal_recursive_matrix(mtxAp)) // fix for serial 20101206
		{
			RSB_DO_ERROR_CUMULATE(errval, rsb__compute_bounded_box(mtxAp));
			goto err;
		}
#pragma omp parallel for schedule(static, 1) reduction(| \
													   : errval) shared(mtxAp) RSB_NTC
		for (smi = 0; smi < mtxAp->all_leaf_matrices_n; ++smi)
		{
			struct rsb_mtx_t *submatrix = mtxAp->all_leaf_matrices[smi].mtxlp;
			RSB_DO_ERROR_CUMULATE(errval, rsb__compute_bounded_box(submatrix));
		}
#pragma omp barrier
	}
	else
	{
#pragma omp parallel for schedule(static, 1) reduction(| \
													   : errval) shared(mtxAp) RSB_NTC
		for (smi = 0; smi < mtxAp->all_leaf_matrices_n; ++smi)
		{
			struct rsb_mtx_t *submatrix = mtxAp->all_leaf_matrices[smi].mtxlp;
			submatrix->bm = submatrix->nr;
			submatrix->bk = submatrix->nc;
			submatrix->broff = submatrix->roff;
			submatrix->bcoff = submatrix->coff;
		}
	}
err:
	RSB_DO_ERR_RETURN(errval)
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
	//printf("Calling rsb_get_num_coo2rec_threads\n");
	// const rsb_thread_t wet = rsb_get_num_coo2rec_threads(); /* want executing threads: */
	const rsb_thread_t wet = nnz; // 1 (global) CUDA thread for each element
	const size_t el_size = rsb__sizeof(typecode);
	rsb_coo_idx_t roff = 0;
	rsb_coo_idx_t coff = 0;
	rsb_nnz_idx_t dnnz = 0;
	//printf("Calling rsb_cuda__estimate_subm_count\n");
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

	//printf("Allocating submatrices\n");
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
	{
		submatricesp[smi] = submatrices + smi;
	}

	//printf("Cleanup\n");
	ect = -dt;
	if ((errval = rsb__do_cleanup_nnz(VA, IA, JA, nnz, roff, coff, m, k, &nnz, typecode, flags)) != RSB_ERR_NO_ERROR)
		goto err;
	ect += (dt = rsb_time());

	est = -dt;
	if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_SORTED_INPUT))
	{
		//printf("Calling rsb__util_sort_row_major_inner\n");
		if ((errval = rsb__util_sort_row_major_inner(VA, IA, JA, nnz, m, k, typecode, flags)) != RSB_ERR_NO_ERROR)
			RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}

	//printf("Adding flags (?)\n");
	RSB_DO_FLAG_ADD(flags, RSB_FLAG_SORTED_INPUT); /* TODO: is this needed ? */
	est += (dt = rsb_time());					   /* with 'sorting' (est) we DO NOT intend also cleanup (in ect) */

	/* we need duplicates removal, and this can only take place after sorting */
	drt = -dt;
	//printf("Calling rsb__weed_out_duplicates\n");
	dnnz = nnz - rsb__weed_out_duplicates(IA, JA, VA, nnz, typecode, flags);
	nnz -= dnnz;
	drt += (dt = rsb_time());
#if defined(RSB__VERBOSE_COO2REC) && RSB__VERBOSE_COO2REC != 0
	RSB_INFO("Duplicates check: %zd - %zd = %zd\n", (size_t)(nnz + dnnz), (size_t)dnnz, (size_t)nnz);
#endif

	/* work vectors allocation */
	/*	IL = rsb__malloc(sizeof(rsb_coo_idx_t)*(m+1)); */
	mat -= dt;
	//printf("Allocating COO vectors\n");
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
	//printf("Calling rsb_cuda_do_compute_vertical_split_parallel\n");
	if ((errval = rsb_cuda_do_compute_vertical_split_parallel(IA, JA, roff, coff, m, k, 0, 0, nnz, IB, NULL, NULL, NULL, NULL, NULL, NULL)) != RSB_ERR_NO_ERROR)
	{
		RSB_PERR_GOTO(err, RSB_ERRM_ES);
	}

#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
	RSB_INFO("beginning (%zd x %zd) @ %p with flags 0x%x (coo:%d, csr:%d), storage: 0x%x, max %zd submatrices\n",
			 (rsb_printf_int_t)submatrices[smi].nr, (rsb_printf_int_t)submatrices[smi].nc, (const void *)&submatrices[smi], submatrices[smi].flags,
			 RSB_DO_FLAG_HAS(submatrices[smi].flags, RSB_FLAG_WANT_COO_STORAGE),
			 RSB_DO_FLAG_HAS(submatrices[smi].flags, RSB_FLAG_WANT_BCSS_STORAGE),
			 submatrices[smi].matrix_storage, (rsb_printf_int_t)tmc);
#endif

	//printf("Calling rsb_do_coo2rec_subdivide_parallel\n");
	errval = rsb_do_coo2rec_subdivide_parallel(VA, IA, JA, m, k, nnz, typecode, pinfop, flags, errvalp, submatricesp, mtxAp, IB, IX, IT, WA, cmc, omc, tmc, RSB_MAX(1, RSB_MIN(wet, nnz)), &cmc);

	sat += (dt = rsb_time());

	//printf("Freeing IX\n");
	RSB_CUDA_CONDITIONAL_FREE(IX);
	if (RSB_SOME_ERROR(errval))
		goto err;

	/*
	RSB_CUDA_CONDITIONAL_FREE(IL);
	RSB_CUDA_CONDITIONAL_FREE(IT);
		*/
	/* WA will is needed for shuffle only (so, after IL,IM deallocation, in a way that total memory need is max(WA,IL)) */
	mat -= dt;
	WA = (rsb_coo_idx_t *)rsb_cuda__malloc(RSB_MAX(sizeof(rsb_coo_idx_t), el_size) * nnz);
	if (!WA)
	{
		errval = RSB_ERR_ENOMEM;
		RSB_PERR_GOTO(err, RSB_ERRM_FAOTAFS);
	}
	mat += (dt = rsb_time());

	eit = -dt;

	for (smi = 0; smi < cmc; ++smi)
	{
		//printf("Calling rsb__is_terminal_recursive_matrix on submatrix %d\n", smi);
		if (rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			++lm;
		}
	}

	/*	qsort(submatricesp+(cmc-lm),(size_t)(lm),sizeof(struct rsb_mtx_t*),&rsb__compar_rcsr_matrix_leftmost_first); */
	qsort(submatricesp, (size_t)(cmc), sizeof(struct rsb_mtx_t *), &rsb__compar_rcsr_matrix_leftmost_first);
	/* TODO: a priority queue would do the job, here */
	for (smi = 0; smi < cmc - lm; ++smi)
	{
		if (rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ANLSMIT);
		}
	}
	for (smi = cmc - lm; smi < cmc; ++smi)
	{
		if (!rsb__is_terminal_recursive_matrix(submatricesp[smi]))
		{
			errval = RSB_ERR_INTERNAL_ERROR;
			RSB_PERR_GOTO(err, RSB_ERRM_ALSMINT);
		}
	}

	//printf("Calling rsb_do_coo2rec_shuffle\n");
	errval = rsb_do_coo2rec_shuffle(VA, IA, JA, m, k, nnz, typecode, pinfop, flags, errvalp, submatricesp, mtxAp, IB, WA, cmc);
	if (RSB_SOME_ERROR(errval))
	{
		RSB_PERR_GOTO(err, RSB_ERRM_SEOWS);
	}

	//printf("Calling rsb__do_set_in_place_submatrices_offsets\n");
	rsb__do_set_in_place_submatrices_offsets(submatrices, cmc, (rsb_char_t *)VA, IA, JA, el_size);

	/*	RSB_INFO("VA:%p, IA:%p, JA:%p\n",VA,IA,JA); */

#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
	RSB_INFO("IA? :%p / %p\n", (const void *)IA,
			 (const void *)(rsb__do_get_first_submatrix(mtxAp)->bpntr -
							(rsb__do_get_first_submatrix(mtxAp)->nr + 1))
			 /*			rsb__do_get_first_submatrix(mtxAp)->roff-
						 ((submatricesp[0])->nr+1) */
	);
#endif

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
	{
#if RSB_WANT_MORE_PARALLELISM
		RSB_DO_ERROR_CUMULATE(errval, rsb_do_switch_fresh_recursive_matrix_to_halfword_storages_parallel(mtxAp));
#else  /* RSB_WANT_MORE_PARALLELISM */
		//printf("Calling rsb_do_switch_fresh_recursive_matrix_to_halfword_storages\n");
		RSB_DO_ERROR_CUMULATE(errval, rsb_do_switch_fresh_recursive_matrix_to_halfword_storages(mtxAp));
#endif /* RSB_WANT_MORE_PARALLELISM */
	}
	else
	{
#if defined(RSB_C2R_IF_VERBOSE) && RSB_C2R_IF_VERBOSE != 0
		RSB_INFO("no  RSB_FLAG_USE_HALFWORD_INDICES flag\n");
#endif
	}

	//printf("Calling rsb_do_compute_bounded_boxes\n");
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
	//printf("Calling rsb__get_index_storage_amount\n");
	mtxAp->idxsa = rsb__get_index_storage_amount(mtxAp);
#endif

	goto noerr;
err:
	mtxAp = NULL;
	RSB_CUDA_CONDITIONAL_FREE(submatrices);
noerr:
	if (RSB_SOME_ERROR(errval))
		rsb__do_perror(NULL, errval);
	RSB_CUDA_CONDITIONAL_FREE(IB);
	RSB_CUDA_CONDITIONAL_FREE(IT);
	RSB_CUDA_CONDITIONAL_FREE(IX);
	/*	RSB_CUDA_CONDITIONAL_FREE(IL); */
	RSB_CUDA_CONDITIONAL_FREE(WA);
	RSB_CUDA_CONDITIONAL_FREE(submatricesp);
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
			//printf("Calling rsb__util_find_coo_max_index_val\n");
			m = rsb__util_find_coo_max_index_val(IA, nnz) + roff + 1;
		}
		if (k == 0 && JA)
		{
			//printf("Calling rsb__util_find_coo_max_index_val\n");
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
		//printf("Calling rsb__do_detect_and_add_triangular_flags\n");
		RSB_DO_FLAG_ADD(flags, rsb__do_detect_and_add_triangular_flags(IA, JA, nnz, flags));
	}

	if (roff && IA)
	{
		//printf("Calling rsb__util_coo_array_add\n");
		rsb__util_coo_array_add(IA, nnz, roff);
	}
	if (coff && JA)
	{
		//printf("Calling rsb__util_coo_array_add\n");
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
			//printf("Calling rsb__allocate_recursive_sparse_matrix_from_row_major_coo\n");
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
		//printf("Calling rsb_cuda__util_coo_alloc_copy_and_stats\n");
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

	//printf("Calling rsb_cuda__mtx_alloc_inner\n");
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
	//printf("calling rsb_cuda__do_mtx_alloc_from_coo_const\n");
	mtxAp = rsb_cuda__do_mtx_alloc_from_coo_const(VA, IA, JA, nnzA, typecode, nrA, ncA, brA, bcA, flagsA, &errval);
	RSB_INTERFACE_RETURN_MTX_ERRP(mtxAp, errval, errvalp);
}

struct rsb_mtx_t *rsb_cuda__do_mtx_free(struct rsb_mtx_t *mtxAp); // this one is here just to make the rsb_cuda__destroy_inner function compile

void *rsb_cuda__destroy_inner(struct rsb_mtx_t *mtxAp)
{
	rsb_submatrix_idx_t i, j;
	struct rsb_mtx_t *submatrix = NULL;

	if (!mtxAp)
	{
		goto ret;
	}

	if (RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS))
	{
		if (rsb__is_root_matrix(mtxAp))
		{
			/* FIXME: unfinished, temporary */
			/* this is a trick: fitting the whole recursive matrix in three arrays */
			void *IA = NULL, *JA = NULL, *VA = NULL;
			const rsb_bool_t is_bio = true; // rsb__do_is_matrix_binary_loaded(mtxAp); // binary I/O matrix
			struct rsb_mtx_t *fsm = rsb__do_get_first_submatrix(mtxAp);
			rsb_flags_t flags = mtxAp->flags;

			if (!is_bio)
			{
				RSB_CUDA_CONDITIONAL_FREE(mtxAp->all_leaf_matrices)
			}
			JA = fsm->bindx;
			VA = fsm->VA;
			if (!is_bio)
			{
				// IA = fsm->bpntr-(fsm->nr+1);
				IA = fsm->bpntr;
			}
			else
			{
				IA = mtxAp;
			}

			if (!is_bio)
			{
				// extra allocation
				//			RSB_INFO("VA:%p, IA:%p, JA:%p\n",VA,IA,JA);
				RSB_CUDA_CONDITIONAL_FREE(mtxAp)
			}

			if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
			{
				RSB_CUDA_CONDITIONAL_FREE(IA); /* these arrays are allowed to be NULL, as it happens during conversions */
				RSB_CUDA_CONDITIONAL_FREE(JA);
				RSB_CUDA_CONDITIONAL_FREE(VA);
			}
		}
		return NULL; /* no deallocation, in this case */
	}

	if (rsb__is_recursive_matrix(mtxAp->flags))
		RSB_SUBMATRIX_FOREACH(mtxAp, submatrix, i, j)
		{
			if (submatrix)
			{
				rsb_cuda__do_mtx_free(submatrix);
			}
		}

	if (!(mtxAp->flags & RSB_FLAG_EXPERIMENTAL_IN_PLACE_CSR))
		if (!RSB_DO_FLAG_HAS(mtxAp->flags, RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS))
		{
			RSB_CUDA_CONDITIONAL_FREE(mtxAp->VA);
			RSB_CUDA_CONDITIONAL_FREE(mtxAp->bindx);
			RSB_CUDA_CONDITIONAL_FREE(mtxAp->bpntr);
		}

	RSB_CUDA_CONDITIONAL_FREE(mtxAp->indptr);

#if RSB_WANT_BITMAP
	if (mtxAp->options)
		rsb__destroy_options_t(mtxAp->options);
#endif /* RSB_WANT_BITMAP */

	if ((mtxAp->flags & RSB_FLAG_OWN_PARTITIONING_ARRAYS) != 0)
	{
		RSB_CUDA_CONDITIONAL_FREE(mtxAp->rpntr);
		RSB_CUDA_CONDITIONAL_FREE(mtxAp->cpntr);
	}
#if RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS
	RSB_CUDA_CONDITIONAL_FREE(mtxAp->all_leaf_matrices);
#endif					/* RSB_EXPERIMENTAL_SHOULD_TRAVERSE_RECURSIVE_MATRICES_AS_BLOCKS */
	RSB_BZERO_P(mtxAp); /* this enforces correct usage */
ret:
	return NULL;
}

struct rsb_mtx_t *rsb_cuda__do_mtx_free(struct rsb_mtx_t *mtxAp)
{
	rsb_flags_t flags;

	if (!mtxAp)
	{
		goto ret;
	}

	/*if (RSB_MTX_HBDF(mtxAp))
	{
		blas_sparse_matrix bmtxA = RSB_MTX_HBDFH(mtxAp);
		rsb__BLAS_Xusds(bmtxA);
		RSB_CUDA_CONDITIONAL_FREE(mtxAp);
		goto ret;
	}*/

	flags = mtxAp->flags;

	rsb_cuda__destroy_inner(mtxAp);

	if (!RSB_DO_FLAG_HAS(flags, RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS))
	{
		if (mtxAp && RSB_WANT_ZERO_ON_DESTROY)
		{
			RSB_BZERO_P(mtxAp);
		}
		RSB_CUDA_CONDITIONAL_FREE(mtxAp);
	}
	else
	{
		mtxAp = NULL;
	}

	RSB_DEBUG_ASSERT(!mtxAp);
ret:
	return mtxAp;
}

struct rsb_mtx_t *rsb_cuda_mtx_free(struct rsb_mtx_t *mtxAp)
{
	struct rsb_mtx_t *mtxBp = NULL;
	RSB_INTERFACE_PREAMBLE
	mtxBp = rsb_cuda__do_mtx_free(mtxAp);
	return mtxBp;
}