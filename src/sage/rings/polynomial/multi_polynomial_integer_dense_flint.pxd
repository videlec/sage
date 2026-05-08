#
# Declarations for fmpz_mpoly (multivariate polynomials over ZZ in FLINT)
#

# Import FLINT types
from sage.libs.flint.flint cimport slong, fmpz_t

# Ordering types for multivariate polynomials
cdef enum ordering_t:
    ORDERING_LEX = 0          # Lexicographic order
    ORDERING_DEGLEX = 1       # Degree lexicographic order
    ORDERING_DEGREVLEX = 2   # Degree reverse lexicographic order

# Opaque type for fmpz_mpoly
cdef extern from "fmpz_mpoly.h":
    ctypedef struct fmpz_mpoly_struct:
        pass
    ctypedef fmpz_mpoly_struct fmpz_mpoly_t[1]

    # Initialization and memory management
    void fmpz_mpoly_init(fmpz_mpoly_t, slong nvars, ordering_t)
    void fmpz_mpoly_clear(fmpz_mpoly_t)
    void fmpz_mpoly_set(fmpz_mpoly_t res, fmpz_mpoly_t a)
    void fmpz_mpoly_swap(fmpz_mpoly_t a, fmpz_mpoly_t b)

    # Basic arithmetic
    void fmpz_mpoly_add(fmpz_mpoly_t res, fmpz_mpoly_t a, fmpz_mpoly_t b)
    void fmpz_mpoly_sub(fmpz_mpoly_t res, fmpz_mpoly_t a, fmpz_mpoly_t b)
    void fmpz_mpoly_mul(fmpz_mpoly_t res, fmpz_mpoly_t a, fmpz_mpoly_t b)
    void fmpz_mpoly_neg(fmpz_mpoly_t res, fmpz_mpoly_t a)
    void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t res, fmpz_mpoly_t a, fmpz_t c)

    # Generators and constants
    void fmpz_mpoly_gen(fmpz_mpoly_t res, slong var, ordering_t)
    void fmpz_mpoly_one(fmpz_mpoly_t res)
    void fmpz_mpoly_zero(fmpz_mpoly_t res)

    # String representation
    char* fmpz_mpoly_get_str(fmpz_mpoly_t a)
    int fmpz_mpoly_set_str(fmpz_mpoly_t poly, const char *str)

    # Degree and coefficient access
    slong fmpz_mpoly_total_degree(fmpz_mpoly_t a)
    slong fmpz_mpoly_length(fmpz_mpoly_t poly)
    void fmpz_mpoly_get_term_exp(fmpz_mpoly_t poly, slong i, slong *exps)
    void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t coeff, fmpz_mpoly_t poly, slong i)
    void fmpz_mpoly_set_coeff_fmpz(fmpz_mpoly_t poly, fmpz_t c, slong *exps)

    # Evaluation
    void fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t res, fmpz_mpoly_t poly, slong var, fmpz_t val)
    void fmpz_mpoly_evaluate_all_fmpz(fmpz_t res, fmpz_mpoly_t poly, fmpz_t *vals, slong *exps)

    # Comparison
    int fmpz_mpoly_equal(fmpz_mpoly_t a, fmpz_mpoly_t b)
    int fmpz_mpoly_is_zero(fmpz_mpoly_t a)
    int fmpz_mpoly_is_one(fmpz_mpoly_t a)

    # GCD
    void fmpz_mpoly_gcd(fmpz_mpoly_t res, fmpz_mpoly_t a, fmpz_mpoly_t b)
