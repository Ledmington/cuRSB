MACROS=-DRSB_CONST_MAX_SUPPORTED_THREADS=128 -DRSB_WANT_VERBOSE_MESSAGES=0
CFILES=librsb/rsb_rsb.c librsb/rsb_init.c librsb/rsb_internals.c librsb/rsb_do.c librsb/rsb_clone.c   \
		librsb/rsb_spmv.c librsb/rsb_spsv.c librsb/rsb_srtp.c librsb/rsb_coo.c librsb/rsb_mio.c       \
		librsb/rsb_err.c librsb/rsb_spsum.c librsb/rsb_spgemm.c librsb/rsb_libspblas_handle.c         \
		librsb/rsb_get.c librsb/rsb_eps.c librsb/rsb_strmif.c librsb/rsb_sys.c librsb/rsb_tune.c      \
		librsb/rsb_swt.c librsb/rsb_util.c librsb/rsb_bio.c librsb/rsb_perf.c librsb/rsb_user.c       \
		librsb/rsb_test_accuracy.c librsb/rsb_is.c librsb/rsb_rec.c librsb/rsb_srt.c librsb/rsb_set.c \
		librsb/rsb_cpmv.c librsb/rsb_asm.c librsb/rsb_stropts.c librsb/rsb_idx.c librsb/rsb_dump.c    \
		librsb/rsb_prec.c librsb/rsb_csr2coo.c librsb/rsb_krnl.c librsb/rsb_spsum_misc.c              \
		librsb/rsb_rec2coo.c librsb/rsb_render.c librsb/rsb_coo2rec.c librsb/rsb_gen.c                \
		librsb/rsb_blas_stuff.c librsb/rsb_permute.c librsb/rsb_merge.c librsb/rsb_msort_up.c         \
		librsb/rsb_mmio.c librsb/rsb_spgemm_csr.c librsb/rsb_src.c librsb/rsb_garbage.c               \
		librsb/rsb_bench.c librsb/rsb_coo_check.c librsb/rsb_csr.c librsb/rsb_krnl_bcoo_spmv_u.c      \
		librsb/rsb_krnl_bcss*.c
GCCFLAGS=-Xcompiler -fopenmp
LINKING=-lgomp

build:
	nvcc ${MACROS} ${GCCFLAGS} librsb/examples/hello_cuda.cu ${CFILES} librsb/*.cu -o librsb/examples/hello_cuda ${LINKING}
	nvcc ${MACROS} ${GCCFLAGS} librsb/examples/rsb_cuda_bench.cu ${CFILES} librsb/*.cu -o librsb/examples/rsb_cuda_bench ${LINKING}
