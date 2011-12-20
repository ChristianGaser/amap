/* CLAPACK 3.0 BLAS wrapper macros
 * Feb 5, 2000
 */

#ifndef __BLASWRAP_H
#define __BLASWRAP_H

#ifndef NO_BLAS_WRAP

#define dsytrf_ EBTKS_dsytrf
#define dsytrs_ EBTKS_dsytrs
#define xerbla_ EBTKS_xerbla
#define lsame_ EBTKS_lsame
#define ilaenv_ EBTKS_ilaenv
#define dsytf2_ EBTKS_dsytf2
#define dlasyf_ EBTKS_dlasyf
#define ieeeck_ EBTKS_ieeeck
#define s_cmp EBTKS_s_cmp
#define s_copy EBTKS_s_copy

#define dswap_ EBTKS_dswap
#define dscal_ EBTKS_dscal
#define dcopy_ EBTKS_dcopy
#define idamax_ EBTKS_idamax
 
/* BLAS2 routines */
#define dgemv_ EBTKS_dgemv
#define dger_ EBTKS_dger
#define dsyr_ EBTKS_dsyr
 
/* BLAS3 routines */
#define dgemm_ EBTKS_dgemm
#define dsyrk_ EBTKS_dsyrk
#define dsyr2k_ EBTKS_dsyr2k
#define dtrmm_ EBTKS_dtrmm
#define dtrsm_ EBTKS_dtrsm

#endif /* NO_BLAS_WRAP */

#endif /* __BLASWRAP_H */
