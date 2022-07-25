#ifndef RSBCK_HPP_INCLUDED
#define RSBCK_HPP_INCLUDED

#include "rsbpp.hpp"

// <double,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<double,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<float,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_coo_idx_t* IP, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <double,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<double,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<float,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_csr_spmx<std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_coo_idx_t* IP, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);

// <double,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<double,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<float,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<float>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<double>,rsb_coo_idx_t,rsb_half_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_half_idx_t* IA, const rsb_half_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <double,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<double,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const double* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const double* rhs, const rsb_coo_idx_t ldY, double* out, const double* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <float,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<float,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const float* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const float* rhs, const rsb_coo_idx_t ldY, float* out, const float* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<float>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<float>* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<float>* rhs, const rsb_coo_idx_t ldY, std::complex<float>* out, const std::complex<float>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);
// <std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>
extern template rsb_err_t rsbpp_coo_spmx<std::complex<double>,rsb_coo_idx_t,rsb_coo_idx_t>(rsb_flags_t flags, const rsb_coo_idx_t nnz, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const std::complex<double>* VA, const rsb_coo_idx_t* IA, const rsb_coo_idx_t* JA, const rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const std::complex<double>* rhs, const rsb_coo_idx_t ldY, std::complex<double>* out, const std::complex<double>* alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, bool by_rows);

#endif /* RSBCK_HPP_INCLUDED */
