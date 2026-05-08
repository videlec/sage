#
# Multivariate Polynomials over ZZ using FLINT's fmpz_mpoly
#

# Import Python and Sage dependencies
from sage.rings.polynomial.multi_polynomial_element import MPolynomialElement
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.integer_ring import ZZ
from sage.structure.parent import Parent
from sage.structure.element import RingElement
from sage.libs.flint.flint_integer cimport fmpz

# Import Cython declarations
from .multi_polynomial_integer_dense_flint cimport (
    fmpz_mpoly_t, ordering_t, ORDERING_LEX, ORDERING_DEGLEX, ORDERING_DEGREVLEX,
    fmpz_mpoly_init, fmpz_mpoly_clear, fmpz_mpoly_add, fmpz_mpoly_sub,
    fmpz_mpoly_mul, fmpz_mpoly_neg, fmpz_mpoly_gen, fmpz_mpoly_one,
    fmpz_mpoly_zero, fmpz_mpoly_is_zero, fmpz_mpoly_is_one, fmpz_mpoly_get_str,
    fmpz_mpoly_set_str, fmpz_mpoly_total_degree, fmpz_mpoly_length,
    fmpz_mpoly_get_term_exp, fmpz_mpoly_get_term_coeff_fmpz, fmpz_mpoly_set_coeff_fmpz,
    fmpz_mpoly_evaluate_one_fmpz, fmpz_mpoly_scalar_mul_fmpz, fmpz_mpoly_equal,
    fmpz_mpoly_set, fmpz_mpoly_gcd
)

# Map Sage ordering strings to FLINT ordering_t
ORDERING_MAP = {
    'lex': ORDERING_LEX,
    'deglex': ORDERING_DEGLEX,
    'degrevlex': ORDERING_DEGREVLEX,
}

# Parent class for the ring
class MPolynomialRing_integer_dense_flint(MPolynomialRing_base):
    def __init__(self, n, names=None, order='lex'):
        """
        Initialize a multivariate polynomial ring over ZZ using FLINT.

        INPUT:

        - ``n``: number of variables.
        - ``names``: list of variable names (default: ['x0', 'x1', ...]).
        - ``order``: monomial ordering (default: 'lex').
        """
        if names is None:
            names = [f'x{i}' for i in range(n)]
        if order not in ORDERING_MAP:
            raise ValueError(f"Unsupported ordering: {order}. Use one of {list(ORDERING_MAP.keys())}")

        self._n = n
        self._names = names
        self._order = order
        self._flint_order = ORDERING_MAP[order]

        super().__init__(ZZ, n, names, order=order)

    def _repr_(self):
        return f"Multivariate Polynomial Ring in {', '.join(self._names)} over Integer Ring (FLINT backend)"

    def base_ring(self):
        return ZZ

    def characteristic(self):
        return 0

    def is_field(self, proof=True):
        return False

    def is_commutative(self):
        return True

    def gen(self, i):
        """
        Return the i-th generator (x_i).
        """
        cdef fmpz_mpoly_t poly
        fmpz_mpoly_init(poly, self._n, self._flint_order)
        fmpz_mpoly_gen(poly, i, self._flint_order)
        return self.element_class(self, poly)

    def gens(self):
        """
        Return the generators of this ring as a tuple.
        """
        return tuple(self.gen(i) for i in range(self._n))

    def ngens(self):
        return self._n

    def variable_names(self):
        return self._names

    def zero(self):
        """
        Return the zero element of this ring.
        """
        cdef fmpz_mpoly_t poly
        fmpz_mpoly_init(poly, self._n, self._flint_order)
        fmpz_mpoly_zero(poly)
        return self.element_class(self, poly)

    def one(self):
        """
        Return the one element of this ring.
        """
        cdef fmpz_mpoly_t poly
        fmpz_mpoly_init(poly, self._n, self._flint_order)
        fmpz_mpoly_one(poly)
        return self.element_class(self, poly)

    def _element_constructor_(self, x, check=True):
        """
        Construct an element of this ring.
        """
        if isinstance(x, int):
            cdef fmpz_mpoly_t poly
            fmpz_mpoly_init(poly, self._n, self._flint_order)
            if x == 0:
                fmpz_mpoly_zero(poly)
            elif x == 1:
                fmpz_mpoly_one(poly)
            else:
                cdef fmpz_t c
                fmpz_init(c)
                fmpz_set_si(c, x)
                fmpz_mpoly_set_coeff_fmpz(poly, c, NULL)
                fmpz_clear(c)
            return self.element_class(self, poly)

        elif isinstance(x, str):
            cdef fmpz_mpoly_t poly
            fmpz_mpoly_init(poly, self._n, self._flint_order)
            if fmpz_mpoly_set_str(poly, x.encode('utf-8')) != 0:
                raise ValueError(f"Could not parse polynomial from string: {x}")
            return self.element_class(self, poly)

        elif isinstance(x, MPolynomial_integer_dense_flint):
            return x

        else:
            raise TypeError(f"Cannot convert {type(x)} to a multivariate polynomial")

    # Element class (defined below)
    Element = None

# Element class for multivariate polynomials
cdef class MPolynomial_integer_dense_flint(MPolynomialElement):
    def __init__(self, parent, poly):
        """
        Initialize a multivariate polynomial over ZZ using FLINT.
        """
        self._parent = parent
        self._poly = poly

    def __dealloc__(self):
        if self._poly is not None:
            fmpz_mpoly_clear(self._poly)

    @staticmethod
    cdef _new_from_fmpz_mpoly(parent, fmpz_mpoly_t poly):
        """
        Create a new element from an fmpz_mpoly_t.
        """
        cdef MPolynomial_integer_dense_flint P = <MPolynomial_integer_dense_flint>Parent.__new__(MPolynomial_integer_dense_flint)
        P._parent = parent
        P._poly = poly
        return P

    def __repr__(self):
        return f"fmpz_mpoly({fmpz_mpoly_get_str(self._poly).decode('utf-8')})"

    def __str__(self):
        return fmpz_mpoly_get_str(self._poly).decode('utf-8')

    def __hash__(self):
        return hash(str(self))

    # Arithmetic operations
    def __add__(self, other):
        """
        Add two polynomials.
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_add(result, self._poly, other._poly)
        return self._parent.element_class(self._parent, result)

    def __sub__(self, other):
        """
        Subtract two polynomials.
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_sub(result, self._poly, other._poly)
        return self._parent.element_class(self._parent, result)

    def __mul__(self, other):
        """
        Multiply two polynomials or by a scalar.
        """
        if isinstance(other, int):
            cdef fmpz_t c
            cdef fmpz_mpoly_t result
            fmpz_init(c)
            fmpz_set_si(c, other)
            fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
            fmpz_mpoly_scalar_mul_fmpz(result, self._poly, c)
            fmpz_clear(c)
            return self._parent.element_class(self._parent, result)
        elif not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_mul(result, self._poly, other._poly)
        return self._parent.element_class(self._parent, result)

    def __rmul__(self, other):
        """
        Right multiplication by a scalar.
        """
        return self.__mul__(other)

    def __neg__(self):
        """
        Negate the polynomial.
        """
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_neg(result, self._poly)
        return self._parent.element_class(self._parent, result)

    def __pow__(self, exponent):
        """
        Raise the polynomial to a power.
        """
        if not isinstance(exponent, int):
            raise TypeError("Exponent must be an integer")
        if exponent < 0:
            raise ValueError("Negative exponent not supported")
        if exponent == 0:
            return self._parent.one()
        if exponent == 1:
            return self

        cdef fmpz_mpoly_t result, temp
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_init(temp, self._parent._n, self._parent._flint_order)

        fmpz_mpoly_set(result, self._poly)
        for _ in range(exponent - 1):
            fmpz_mpoly_mul(temp, result, self._poly)
            fmpz_mpoly_set(result, temp)

        fmpz_mpoly_clear(temp)
        return self._parent.element_class(self._parent, result)

    # Comparison
    def __eq__(self, other):
        """
        Check equality with another polynomial.
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            return False
        return bool(fmpz_mpoly_equal(self._poly, other._poly))

    def __ne__(self, other):
        """
        Check inequality with another polynomial.
        """
        return not self.__eq__(other)

    # Properties
    def is_zero(self):
        """
        Check if the polynomial is zero.
        """
        return bool(fmpz_mpoly_is_zero(self._poly))

    def is_one(self):
        """
        Check if the polynomial is one.
        """
        return bool(fmpz_mpoly_is_one(self._poly))

    def total_degree(self):
        """
        Return the total degree of the polynomial.
        """
        return fmpz_mpoly_total_degree(self._poly)

    def monomials(self):
        """
        Return the list of monomials of this polynomial.
        """
        cdef slong i, length = fmpz_mpoly_length(self._poly)
        cdef slong *exps = <slong*>malloc(self._parent._n * sizeof(slong))
        monomials = []
        for i in range(length):
            fmpz_mpoly_get_term_exp(self._poly, i, exps)
            monomial = []
            for j in range(self._parent._n):
                monomial.append(exps[j])
            monomials.append(tuple(monomial))
        free(exps)
        return monomials

    def coefficients(self):
        """
        Return the list of coefficients of this polynomial.
        """
        cdef slong i, length = fmpz_mpoly_length(self._poly)
        cdef fmpz_t coeff
        cdef slong *exps = <slong*>malloc(self._parent._n * sizeof(slong))
        coefficients = []
        fmpz_init(coeff)
        for i in range(length):
            fmpz_mpoly_get_term_coeff_fmpz(coeff, self._poly, i)
            coefficients.append(int(fmpz_get_si(coeff)))
        fmpz_clear(coeff)
        free(exps)
        return coefficients

    def subs(self, **kwargs):
        """
        Substitute variables with values.
        """
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_set(result, self._poly)

        for var, value in kwargs.items():
            if not isinstance(value, int):
                raise TypeError(f"Substitution value for {var} must be an integer")
            try:
                idx = self._parent._names.index(var)
            except ValueError:
                raise ValueError(f"Variable {var} not in the ring")

            cdef fmpz_t val
            fmpz_init(val)
            fmpz_set_si(val, value)
            fmpz_mpoly_evaluate_one_fmpz(result, result, idx, val)
            fmpz_clear(val)

        return self._parent.element_class(self._parent, result)

    def gcd(self, other):
        """
        Compute the GCD of this polynomial and another.
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_gcd(result, self._poly, other._poly)
        return self._parent.element_class(self._parent, result)

# Assign the Element class to the parent
MPolynomialRing_integer_dense_flint.Element = MPolynomial_integer_dense_flint
