#include <rsb.h>         /* for rsb_lib_init */
#include <blas_sparse.h> /* Sparse BLAS on the top of librsb */

#define RSB_CHECK_ERROR(X) if(X != RSB_ERR_NO_ERROR){return -1;}

int main(const int argc, char * const argv[]) {
	blas_sparse_matrix matrix = blas_invalid_handle; /* handle for A */
	const int nnzA = 4;	/* number of nonzeroes of matrix A */
	const int nrA = 3, ncA = 3;	/* number of A's rows and columns */
	const int incX = 1, incB = 1; /* spacing of X, B entries */
	int IA[] = { 0, 1, 2, 2 }; /* row indices*/
	int JA[] = { 0, 1, 0, 2 }; /* column indices*/
	double VA[] = { 11.0, 22.0, 13.0, 33.0  }; /* coefficients */
	double X[] = {  0.0,  0.0,  0.0 };
	double B[] = { -1.0, -2.0, -2.0 };
	rsb_err_t errval = RSB_ERR_NO_ERROR;

    // Initialize librsb
    RSB_CHECK_ERROR(rsb_lib_init(RSB_NULL_INIT_OPTIONS))

	matrix = BLAS_duscr_begin(nrA, ncA);                             /* begin matrix creation */
	if (matrix == blas_invalid_handle) return -1;                    /* check */
	if ( BLAS_ussp(matrix,blas_lower_symmetric) != 0 ) return -1;  /* set symmetry property */
	if ( BLAS_duscr_insert_entries(matrix, nnzA, VA, IA, JA) != 0 ) return -1; /* insert some nonzeroes */
	if ( BLAS_duscr_end(matrix) == blas_invalid_handle ) return -1; /* end creation */

	if ( BLAS_dusmv(blas_no_trans, -1, matrix, B, incB, X, incX)) {
        return -1;    /* X := X + (-1) * A * B  */
    }

	if ( BLAS_usds(matrix)) {
        return -1;   // destroy matrix
    }
	
    errval = rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
    if (errval != RSB_ERR_NO_ERROR) {
        return -1;
    }

	return 0;
}