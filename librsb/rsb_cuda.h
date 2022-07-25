#ifndef RSB_CUDA_H_INCLUDED
#define RSB_CUDA_H_INCLUDED

rsb_err_t rsb_cuda_spmv(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, const void * Xp, rsb_coo_idx_t incX, const void * betap, void * Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_cuda_spmm(rsb_trans_t transA, const void * alphap, const struct rsb_mtx_t * mtxAp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * Bp, rsb_nnz_idx_t ldB, const void * betap, void * Cp, rsb_nnz_idx_t ldC);
rsb_err_t rsb_cuda_spsv(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, const void * Xp, rsb_coo_idx_t incX, void * Yp, rsb_coo_idx_t incY);
rsb_err_t rsb_cuda_spsm(rsb_trans_t transT, const void * alphap, const struct rsb_mtx_t * mtxTp, rsb_coo_idx_t nrhs, rsb_flags_t order, const void * betap, const void * Bp, rsb_nnz_idx_t ldB, void * Cp, rsb_nnz_idx_t ldC);

#endif // RSB_CUDA_H_INCLUDED