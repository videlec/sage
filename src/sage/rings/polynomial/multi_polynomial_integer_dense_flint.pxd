#
# Declarations for multivariate polynomials over ZZ using FLINT's fmpz_mpoly
#

# Import existing FLINT types and functions from SageMath
from sage.libs.flint.types cimport fmpz_t, fmpz_mpoly_t, ordering_t, ORDERING_LEX, ORDERING_DEGLEX, ORDERING_DEGREVLEX
from sage.libs.flint.fmpz_mpoly cimport (
    fmpz_mpoly_init, fmpz_mpoly_clear, fmpz_mpoly_set, fmpz_mpoly_swap,
    fmpz_mpoly_add, fmpz_mpoly_sub, fmpz_mpoly_mul, fmpz_mpoly_neg,
    fmpz_mpoly_gen, fmpz_mpoly_one, fmpz_mpoly_zero,
    fmpz_mpoly_get_str, fmpz_mpoly_set_str,
    fmpz_mpoly_total_degree, fmpz_mpoly_length,
    fmpz_mpoly_get_term_exp, fmpz_mpoly_get_term_coeff_fmpz, fmpz_mpoly_set_coeff_fmpz,
    fmpz_mpoly_evaluate_one_fmpz, fmpz_mpoly_scalar_mul_fmpz,
    fmpz_mpoly_equal, fmpz_mpoly_is_zero, fmpz_mpoly_is_one,
    fmpz_mpoly_gcd
)
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_set_si, fmpz_set_mpz, fmpz_get_mpz

# Declare the parent class for the ring
cdef class MPolynomialRing_integer_dense_flint

# Declare the element class
cdef class MPolynomial_integer_dense_flint
