dnl
dnl
dnl	@author: Michele Martone
dnl
/* @cond INNERDOC */
include(`rsb_krnl_macros.m4')dnl
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains functions for sparse recursive multicore matrix vector multiplication.
 */
/*
 * FIXME: many beta-related ops are NOT parallel, and this is BAD.
 *
 * */
RSB_M4_HEADER_MESSAGE()dnl


dnl
ifdef(`ONLY_WANT_HEADERS',`
#ifndef RSB_SPMV_H_INCLUDED
#define RSB_SPMV_H_INCLUDED
')
dnl
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
dnl
#include "rsb_internals.h"
#include "rsb_lock.h"
dnl 
ifdef(`ONLY_WANT_HEADERS',`dnl
dnl
#define RSB_ENABLE_INNER_NRHS_SPMV 1
dnl
#if RSB_ENABLE_INNER_NRHS_SPMV
#define RSB_INNER_NRHS_SPMV_ARGS	,const rsb_int_t nrhs, /*const size_t outtot, const size_t rhstot,*/ const size_t outnri, const size_t rhsnri
#define RSB_INNER_NRHS_SPMV_ARGS_IDS	,nrhs/*,outtot,rhstot*/,outnri,rhsnri
#define RSB_INNER_NRHS_SPMV_YSCA_IDS	,nrhs,outnri
#define RSB_OUTER_NRHS_SPMV_ARGS	,const rsb_int_t nrhs, const size_t outnri, const size_t rhsnri
#define RSB_OUTER_NRHS_SPMV_ARGS_IDS	,nrhs,outnri,rhsnri
#else /* RSB_ENABLE_INNER_NRHS_SPMV */
#define RSB_INNER_NRHS_SPMV_ARGS	
#define RSB_INNER_NRHS_SPMV_YSCA_IDS		/* */
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
#define RSB_DEFAULT_INNER_NRHS_SPMV_ARGS	,1,/*0,0,*/0,0
#define RSB_DEFAULT_OUTER_NRHS_SPMV_ARGS	,1,0,0
dnl
',`dnl
RSB_INTERNALS_COMMON_HEAD_DECLS
dnl
dnl
/* FIXME: to move these macros to one header in order to avoid any identifier clash */
/* #define RSB_CBLAS_X_SCAL_SPMV rsb__cblas_Xscal */
#define RSB__FOREACH_NRHS(NRHSI, NRHS) for (NRHSI=0;NRHSI<NRHS;++NRHSI)
#define RSB_CBLAS_X_SCAL_SPMV(TYPECODE,N,ALPHAP,A,STRIDE) rsb__cblas_Xscal_parallel((TYPECODE),(N),(ALPHAP),(A),(STRIDE))
#if RSB_ENABLE_INNER_NRHS_SPMV
#define RSB_CBLAS_X_SCAL_SPMM(TYPECODE,N,ALPHAP,A,STRIDE) 					\
{	/* FIXME: this interacts with RSB_INNER_NRHS_SPMV_ARGS */					\
	rsb_int_t nrhsi = 0;										\
	RSB__FOREACH_NRHS(nrhsi,nrhs)									\
	{												\
		RSB_CBLAS_X_SCAL_SPMV(TYPECODE, N, ALPHAP, RSB_TYPED_OFF_PTR(TYPECODE,A,nrhsi*(outnri)), STRIDE); 	\
	}												\
}										/* RSB_CBLAS_X_SCAL_SPMM */
#else  /* RSB_ENABLE_INNER_NRHS_SPMV */
#define RSB_CBLAS_X_SCAL_SPMM(TYPECODE,N,ALPHAP,A,STRIDE) RSB_CBLAS_X_SCAL_SPMV(TYPECODE,N,ALPHAP,A,STRIDE) 
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
')dnl

dnl
rsb_err_t rsb_do_spmv_non_recursive(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_err_t errval = RSB_ERR_NO_ERROR;
	rsb_bool_t nostride = ( incx == 1 && incy == 1 )?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	rsb_bool_t should_scale_y = ( betap && !RSB_IS_ELEMENT_ONE( betap,mtxAp->typecode))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	rsb_bool_t use_alpha_one = (!alphap || RSB_IS_ELEMENT_ONE(alphap,mtxAp->typecode))?RSB_BOOL_TRUE:RSB_BOOL_FALSE;
	rsb_bool_t use_y_zeroing_kernel = (should_scale_y && RSB_IS_ELEMENT_ZERO(betap,mtxAp->typecode) && nostride && use_alpha_one)?RSB_BOOL_TRUE:RSB_BOOL_FALSE;

	/*
		FIXME : should handle beta in a more specialized fashion.
			should also handle more specialized alphap cases.
	*/
#if RSB_ENABLE_INNER_NRHS_SPMV
	/*const size_t outtot=0,rhstot=0;
	const size_t outnri=0,rhsnri=0;
	const rsb_int_t nrhs=1;*/
 	/* FIXME: the above should be specified from argument */
	const size_t lenx=(mtxAp->el_size*rhsnri);
	const size_t leny=(mtxAp->el_size*outnri);
	rsb_int_t nrhsi=0;
	for(nrhsi=0;nrhsi<nrhs;++nrhsi)
	{
		void      *out=((      rsb_byte_t*)y)+(leny*nrhsi);
		const void*rhs=((const rsb_byte_t*)x)+(lenx*nrhsi);
#else /* RSB_ENABLE_INNER_NRHS_SPMV */
		void      *out=((      rsb_byte_t*)y);
		const void*rhs=((const rsb_byte_t*)x);
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */

	if(should_scale_y && !use_y_zeroing_kernel)
		RSB_CBLAS_X_SCAL_SPMV(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transA),betap,out,incy);
	/* no beta specified counts as beta=1, and so no scaling is needed */

dnl
dnl		FIXME : yes. using RSB_M4_ARGS_TO_ACTUAL_ARGS two times. we were forced to do so, probably due to a bug in RSB_M4_ARGS_TO_ACTUAL_ARGS.
dnl
	if(use_alpha_one)
	{
		/* no alpha specified counts as alpha=1 */
		if(nostride)
		{
			if(use_y_zeroing_kernel)
				/* y <- a * x  */
				RSB_DO_ERROR_CUMULATE(errval,RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(`spmv_uauz')(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(`spmv_uauz'))))));
			else
				/* y <- y + a * x  */
				RSB_DO_ERROR_CUMULATE(errval,RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(`spmv_uaua')(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(`spmv_uaua'))))));
		}
		else
			/* y <- a * x  , with stride */
			RSB_DO_ERROR_CUMULATE(errval,RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(`spmv_sasa')(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(`spmv_sasa'))))));
	}
	else
	{
		if(nostride)
		{
			/* y <- - a * x  */
			if(RSB_IS_ELEMENT_MINUS_ONE(alphap,mtxAp->typecode))
				RSB_DO_ERROR_CUMULATE(errval,RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(`spmv_unua')(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(`spmv_unua'))))));
			/* y <- alpha * a * x  */
			else
				RSB_DO_ERROR_CUMULATE(errval,RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(`spmv_uxua')(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(`spmv_uxua'))))));
		}
		else
			/* y <- alpha * a * x  , with stride */
			RSB_DO_ERROR_CUMULATE(errval,RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_IDENTIFIER(`spmv_sxsa')(RSB_M4_ARGS_TO_ACTUAL_ARGS(RSB_M4_ARGS_TO_ACTUAL_ARGS((RSB_M4_DIRECT_KERNEL_DISPATCH_FUNCTION_ARGS(`spmv_sxsa'))))));
	}
dnl
dnl	FIXME : we deliberately ignore other useful kernels we have :
dnl	extra_blas_matrix_ops=spmv_sxsx,spmv_uauz,spmv_uxux
dnl
#if RSB_ENABLE_INNER_NRHS_SPMV
	}
#endif /* RSB_ENABLE_INNER_NRHS_SPMV */
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`',`dnl
#if 0
#define RSB_SPMV_VS_DECL	int*mivr=NULL,mivi=0;
#define RSB_SPMV_VS_ALLOC(MTXAP,ERRVAL,ERRL)	op_flags|=RSB_OP_FLAG_WANT_TRACE_PLOT;if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT){ mivr=rsb__calloc(sizeof(int)*((MTXAP)->all_leaf_matrices_n));if(!mivr){ERRVAL=RSB_ERR_ENOMEM;goto ERRL;} }
#define RSB_SPMV_VS_MARK(MI)		if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT){ mivr[mivi++]=(MI);}
#define RSB_SPMV_VS_DUMP(MTXAP)		if(op_flags & RSB_OP_FLAG_WANT_TRACE_PLOT){ /* for(mivi=0;mivi<((MTXAP)->all_leaf_matrices_n);++mivi) printf("%d ",mivr[mivi]);printf("\n"); */ rsb__dump_postscript_recursion_from_mtx_t(NULL,"spmv-dump.eps",(MTXAP),1,1,512,512,RSB_FLAG_NOFLAGS,0,1,0,mivr); }
#define RSB_SPMV_VS_DEALLOC		RSB_CONDITIONAL_FREE(mivr);
#else /* 0 */
#define RSB_SPMV_VS_DECL
#define RSB_SPMV_VS_ALLOC(MTXAP,ERRVAL,ERRL)
#define RSB_SPMV_VS_MARK(MI)
#define RSB_SPMV_VS_DUMP(MTXAP)
#define RSB_SPMV_VS_DEALLOC
#endif /* 0 */
')dnl

rsb_err_t rsb_do_spmv_recursive_parallel(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA, enum rsb_op_flags_t op_flags RSB_INNER_NRHS_SPMV_ARGS)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	  	\ingroup gr_internals
	*/

	rsb_err_t errval = RSB_ERR_NO_ERROR;
	RSB_SPMV_VS_DECL
#if RSB_WANT_OMP_RECURSIVE_KERNELS

	const struct rsb_translated_matrix_t * all_leaf_matrices=NULL;
	struct rsb_spmv_lock_struct_t lock;
	rsb_submatrix_idx_t all_leaf_matrices_n=0;

dnl	if(alphap || betap || incx>1 || incy>1 || transA != RSB_TRANSPOSITION_N)	/* FIXME */
dnl	{errval = RSB_ERR_UNIMPLEMENTED_YET;goto err;}

	if(!rsb__is_recursive_matrix(mtxAp->flags))
		return rsb_do_spmv_non_recursive(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS);

	all_leaf_matrices  =mtxAp->all_leaf_matrices;
	all_leaf_matrices_n=mtxAp->all_leaf_matrices_n;

	if(!all_leaf_matrices || all_leaf_matrices_n<1)
	{errval = RSB_ERR_ENOMEM;goto err;}

	errval = rsb_do_spmv_lock_init(&lock,rsb_global_session_handle.rsb_want_threads,all_leaf_matrices_n,mtxAp,op_flags,transA,y,incy);
	if(RSB_SOME_ERROR(errval))
		goto err;

#if 0
	if(betap && !RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
	{
	#pragma omp parallel shared(y,mtxAp,rsb_global_session_handle)  RSB_NTC 
	{
		rsb_nnz_idx_t tdim = rsb__do_get_rows_of(mtxAp,transA),dim,chunk;
		rsb_char_t * yy=y;
		rsb_thr_t th_id = omp_get_thread_num();
		if(th_id >= rsb_global_session_handle.rsb_want_threads)
			goto scaled;
		chunk=tdim/rsb_global_session_handle.rsb_want_threads;
		yy+=th_id*chunk*incy;
		if(th_id == rsb_global_session_handle.rsb_want_threads-1)
			dim=tdim-th_id*chunk;
		else		
			dim=chunk;

		if(RSB_IS_ELEMENT_ZERO(betap,mtxAp->typecode))
			RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,dim,NULL,yy,incy);
		else
			RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,dim,betap,yy,incy);
		scaled:
		RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	}
	}
	#pragma omp barrier
#else /* 0 */
	/* TODO: make the following parallel */
	if(betap && !RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
		RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transA),betap,y,incy);
#endif /* 0 */
	RSB_SPMV_VS_ALLOC(mtxAp,errval,err)
	#pragma omp parallel reduction(|:errval) shared(lock,all_leaf_matrices,mtxAp)  RSB_NTC 
{
	rsb_thr_t th_id = omp_get_thread_num();
	rsb_submatrix_idx_t n=0;
	rsb_submatrix_idx_t dm=0;

	if(th_id >= rsb_global_session_handle.rsb_want_threads)
		goto skip;

	if(th_id>=all_leaf_matrices_n)
		goto skip;

again:
	for(n=0;RSB_LIKELY(n<all_leaf_matrices_n);++n)
dnl	//if(!RSB_BITMAP_GET(lock.bmap,1,lock.subms,0,n))
#if RSB_WANT_SM_TO_THREAD_MOD_MAPPING && !RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV 
	if( (n % RSB_MIN(all_leaf_matrices_n,rsb_global_session_handle.rsb_want_threads)) == th_id )
#endif /* RSB_WANT_SM_TO_THREAD_MOD_MAPPING */
	{
		const struct rsb_mtx_t *submatrix=all_leaf_matrices[n].mtxlp;
		char *ov=y;
		rsb_bool_t gomv = RSB_BOOL_FALSE;
		rsb_coo_idx_t oincy=incy;
#if RSB_WANT_SPMV_WITH_REDUCE
		rsb_coo_idx_t rh,r0;	/* new */
#endif /* RSB_WANT_SPMV_WITH_REDUCE */
		#pragma omp critical (rsb_spmv_crs)
#if RSB_WANT_BOUNDED_BOXES_SPMV
		{ gomv=(rsb_do_spmv_lock_get(&lock,th_id,submatrix->broff,submatrix->bm,submatrix->bcoff,submatrix->bk,n,transA,&ov,&oincy)==RSB_BOOL_TRUE); if(gomv==RSB_BOOL_TRUE){RSB_SPMV_VS_MARK(n);} }
#else /* RSB_WANT_BOUNDED_BOXES_SPMV */
		{ gomv=(rsb_do_spmv_lock_get(&lock,th_id,submatrix->roff,submatrix->nr,submatrix->coff,submatrix->nc,n,transA,&ov,&oincy)==RSB_BOOL_TRUE); if(gomv==RSB_BOOL_TRUE){RSB_SPMV_VS_MARK(n);} }
#endif /* RSB_WANT_BOUNDED_BOXES_SPMV */
		if(gomv == RSB_BOOL_TRUE)
		{
			const char * offx=NULL; char *offy=NULL;
			const size_t scoff=submatrix->coff-mtxAp->coff;
			const size_t sroff=submatrix->roff-mtxAp->roff;
dnl			offy=((char*)y)+(mtxAp->el_size*sroff)*incy,offx=((const char*)x)+(mtxAp->el_size*scoff)*incx;
			offy=((char*)ov)+(mtxAp->el_size*sroff)*oincy,offx=((const char*)x)+(mtxAp->el_size*scoff)*incx;

			/* FIXME */
			RSB_ASSERT(scoff>=0);
			RSB_ASSERT(sroff>=0);
dnl #if 1
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_non_recursive(submatrix,offx,offy,alphap,NULL,incx,oincy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
dnl #else
dnl			RSB_DO_ERROR_CUMULATE(errval,rsb__do_spmv_uaua(submatrix,offx,offy,transA));
dnl #endif
                       	#pragma omp critical (rsb_spmv_crs)
			{rsb_do_spmv_lock_release(&lock,th_id,ov);RSB_DO_SPMV_LOCK_DM_INC(lock);}
		}
#if RSB_WANT_SPMV_WITH_REDUCE
		if(gomv == RSB_BOOL_ALMOST_TRUE)
		{
                       	#pragma omp critical (rsb_spmv_crs)
			{rsb__do_pick_candidate_interval_for_reduce(&lock,th_id,&ov,&r0,&rh);}

			if(ov && ov!=y)
			{
				rsb__vectors_left_sum_reduce_and_zero(y,ov,mtxAp->typecode,rh,oincy,r0);
                       		#pragma omp critical (rsb_spmv_crs)
                       		{ rsb__do_release_candidate_interval_for_reduce(&lock,th_id,ov,r0,rh);}
			}
		}
#endif /* RSB_WANT_SPMV_WITH_REDUCE*/
	}
		#pragma omp critical (rsb_spmv_crs)
		{ dm = RSB_DO_SPMV_LOCK_DM(lock); }
		if(dm<all_leaf_matrices_n
#if RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV
			&& ((all_leaf_matrices_n-dm)>th_id)
#endif	/* RSB_WANT_EARLY_PARALLEL_REGION_JUMPOUT_SPMV */
		)goto again;
skip:
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	/* done */
}
	RSB_SPMV_VS_DUMP(mtxAp)
err:
	RSB_SPMV_VS_DEALLOC
#if   !defined(__xlC__)
	RSB_NULL_STATEMENT_FOR_COMPILER_HAPPINESS;
	/* FIXME: xlc does not allow this, but we have experienced problems, without */
	#pragma omp barrier
#endif /* __xlC__ */
	RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_lock_free(&lock));
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = RSB_ERR_UNIMPLEMENTED_YET;
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl
dnl

dnl
rsb_err_t rsb_do_spmv_recursive_serial(const struct rsb_mtx_t * mtxAp, const void * x, void * y, const void *alphap, const void * betap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, rsb_trans_t transA RSB_INNER_NRHS_SPMV_ARGS)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	  	\ingroup gr_internals
		This function does not offer result vector accumulation in case of diagonal implicit matrices.
	*/
	struct rsb_mtx_t * submatrix=NULL;
	rsb_submatrix_idx_t i,j;
	rsb_err_t errval = RSB_ERR_NO_ERROR;

	if( rsb__is_recursive_matrix(mtxAp->flags))
	{
		void*offy=NULL;
		const void *offx=NULL;

		if(betap && !RSB_IS_ELEMENT_ONE(betap,mtxAp->typecode))
		{
			/* should scale the output vector */
			RSB_CBLAS_X_SCAL_SPMM(mtxAp->typecode,rsb__do_get_rows_of(mtxAp,transA),betap,y,incy);
		}

		RSB_SUBMATRIX_FOREACH(mtxAp,submatrix,i,j)
		if(submatrix)
		{
			size_t scoff=submatrix->coff-mtxAp->coff;
			size_t sroff=submatrix->roff-mtxAp->roff;

			/* FIXME */
			RSB_ASSERT(scoff>=0);
			RSB_ASSERT(sroff>=0);

			offy=((char*)y)+(mtxAp->el_size*sroff)*incy,offx=((const char*)x)+(mtxAp->el_size*scoff)*incx;
			RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_recursive_serial(submatrix,offx,offy,alphap,NULL,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
		}
	}
	else
	{
		RSB_DO_ERROR_CUMULATE(errval,rsb_do_spmv_non_recursive(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS));
	}
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl
dnl

dnl
rsb_err_t rsb__do_spmv_general(rsb_trans_t transA, const void *alphap, const struct rsb_mtx_t * mtxAp, const void * x, rsb_coo_idx_t incx, const void * betap, void * y, rsb_coo_idx_t incy, enum rsb_op_flags_t op_flags RSB_OUTER_NRHS_SPMV_ARGS)
ifdef(`ONLY_WANT_HEADERS',`;',`dnl
{
	/**
	  	\ingroup gr_internals
	*/
	rsb_err_t errval = RSB_ERR_BADARGS;

#if RSB_ALLOW_ZERO_DIM
	if(RSB_ANY_MTX_DIM_ZERO(mtxAp))
	{
		errval = RSB_ERR_NO_ERROR;
		goto err; /* FIXME: skipping further checks */
	}
#endif
	if(x==y)
		goto err;

	if(incx<1 || incy<1)
		goto err;

/*	we tolerate NULL alhap and betap */
#if 0
	if(!alphap || !betap)
		goto err;
#endif /* 0 */

	if(!mtxAp || !x || !y || transA == RSB_INVALID_FLAGS)
		goto err;

#if 0
	errval = rsb_do_spmv_recursive_serial(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS);
	goto done;
#endif /* 0 */

#if RSB_WANT_OMP_RECURSIVE_KERNELS
	if(RSB_UNLIKELY(op_flags == RSB_OP_FLAG_WANT_SERIAL))
		errval = rsb_do_spmv_recursive_serial(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS);
	else
	{
		RSB_NUM_THREADS_DECL
		RSB_NUM_THREADS_PUSH
		errval = rsb_do_spmv_recursive_parallel(mtxAp,x,y,alphap,betap,incx,incy,transA,op_flags RSB_OUTER_NRHS_SPMV_ARGS_IDS	);
		RSB_NUM_THREADS_POP
	}
	/* the RSB_OP_FLAG_FAKE_LOCK case is handled by rsb_do_spmv_recursive_parallel */
#else /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	errval = rsb_do_spmv_recursive_serial(mtxAp,x,y,alphap,betap,incx,incy,transA RSB_INNER_NRHS_SPMV_ARGS_IDS);
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
	goto done;
done:
	if(!RSB_UNLIKELY(op_flags&RSB_OP_FLAG_DIAGONAL_OVERRIDE_EXPLICIT))// NEW: fix for odd spsv/diagonal implicit/no-parallel cases
	if(RSB_DO_FLAG_HAS(mtxAp->flags,RSB_FLAG_UNIT_DIAG_IMPLICIT))
	{
		rsb_int_t nrhsi,di;
		const rsb_coo_idx_t ndy = RSB_MIN(mtxAp->nr,mtxAp->nc);
		const int row_major = ( nrhs > 1 && (incx >= nrhs || incy >= nrhs ) );

		if( row_major )
		for (di=0;di<ndy;++di)
			rsb__cblas_Xaxpy(mtxAp->typecode,nrhs,alphap
				,RSB_TYPED_OFF_PTR(mtxAp->typecode,x,incx*di)
				,1
				,RSB_TYPED_OFF_PTR(mtxAp->typecode,y,incy*di)
				,1);
		else
		for (nrhsi=0;nrhsi<nrhs;++nrhsi)
			rsb__BLAS_Xaxpy_parallel(RSB_MIN(mtxAp->nr,mtxAp->nc),alphap
				,RSB_TYPED_OFF_PTR(mtxAp->typecode,y,incy*nrhsi*(rsb__do_get_rows_of(mtxAp,transA)))
				,incy
				,RSB_TYPED_OFF_PTR(mtxAp->typecode,x,incx*nrhsi*(rsb__do_get_columns_of(mtxAp,transA)))
				,incx,mtxAp->typecode);
	}
err:
	RSB_DO_ERR_RETURN(errval)
}
')dnl
dnl
dnl


dnl
#ifdef __cplusplus
}
#endif  /* __cplusplus */
dnl
dnl
ifdef(`ONLY_WANT_HEADERS',`
#endif /* RSB_SPMV_H_INCLUDED */
')
dnl
/* @endcond */
dnl
