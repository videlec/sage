r"""
Multivariate Polynomials over the integers using FLINT's fmpz_mpoly
"""
# ****************************************************************************
#       Copyright (C) 2026 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.multi_polynomial_element cimport MPolynomialElement
from sage.rings.polynomial.multi_polynomial_ring cimport MPolynomialRing_base
from sage.rings.integer_ring cimport ZZ, Integer
from sage.structure.parent cimport Parent
from sage.structure.element cimport RingElement

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

ORDERING_MAP = {
    'lex': ORDERING_LEX,
    'deglex': ORDERING_DEGLEX,
    'degrevlex': ORDERING_DEGREVLEX,
}

class MPolynomialRing_integer_dense_flint(MPolynomialRing_base):
    """
    A multivariate polynomial ring over the integers using FLINT's ``fmpz_mpoly``.

    This class provides a SageMath interface to FLINT's multivariate polynomial
    functionality over the integers.

    EXAMPLES::

        sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
        sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
        sage: R
        Multivariate Polynomial Ring in x, y over Integer Ring (FLINT backend)
        sage: x, y = R.gens()
        sage: x + y
        x + y
    """
    Element = MPolynomial_integer_dense_flint

    def __init__(self, n, names=None, order='lex'):
        """
        Initialize a multivariate polynomial ring over ZZ using FLINT.

        INPUT:

        - ``n``: number of variables.
        - ``names``: list of variable names (default: ['x0', 'x1', ...]).
        - ``order``: monomial ordering (default: 'lex').

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'], order='deglex')
            sage: R
            Multivariate Polynomial Ring in x, y over Integer Ring (FLINT backend)
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
        """
        Return a string representation of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: repr(R)
            'Multivariate Polynomial Ring in x, y over Integer Ring (FLINT backend)'
        """
        return f"Multivariate Polynomial Ring in {', '.join(self._names)} over Integer Ring (FLINT backend)"

    def base_ring(self):
        """
        Return the base ring of this polynomial ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.base_ring()
            Integer Ring
        """
        return ZZ

    def characteristic(self):
        """
        Return the characteristic of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.characteristic()
            0
        """
        return 0

    def is_field(self, proof=True):
        """
        Return whether this ring is a field.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.is_field()
            False
        """
        return False

    def is_commutative(self):
        """
        Return whether this ring is commutative.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.is_commutative()
            True
        """
        return True

    cdef inline _new_element(self, fmpz_mpoly_t poly):
        """
        Create a new element from an fmpz_mpoly_t (inline function).
        """
        cdef MPolynomial_integer_dense_flint P = <MPolynomial_integer_dense_flint>Parent.__new__(MPolynomial_integer_dense_flint)
        P._parent = self
        P._poly = poly
        return P

    def gen(self, i):
        """
        Return the i-th generator of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.gen(0)
            x
            sage: R.gen(1)
            y
        """
        cdef fmpz_mpoly_t poly
        fmpz_mpoly_init(poly, self._n, self._flint_order)
        fmpz_mpoly_gen(poly, i, self._flint_order)
        return self._new_element(poly)

    def gens(self):
        """
        Return the generators of this ring as a tuple.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.gens()
            (x, y)
        """
        return tuple(self.gen(i) for i in range(self._n))

    def ngens(self):
        """
        Return the number of generators of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.ngens()
            2
        """
        return self._n

    def variable_names(self):
        """
        Return the names of the variables of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.variable_names()
            ['x', 'y']
        """
        return self._names

    def zero(self):
        """
        Return the zero element of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.zero()
            0
        """
        cdef fmpz_mpoly_t poly
        fmpz_mpoly_init(poly, self._n, self._flint_order)
        fmpz_mpoly_zero(poly)
        return self._new_element(poly)

    def one(self):
        """
        Return the one element of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.one()
            1
        """
        cdef fmpz_mpoly_t poly
        fmpz_mpoly_init(poly, self._n, self._flint_order)
        fmpz_mpoly_one(poly)
        return self._new_element(poly)

    def _element_constructor_(self, x, check=True):
        """
        Construct an element of this ring.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R(1)
            1
            sage: R(0)
            0
            sage: R(Integer(5))
            5
        """
        if isinstance(x, (int, Integer)):
            cdef fmpz_mpoly_t poly
            fmpz_mpoly_init(poly, self._n, self._flint_order)
            if x == 0:
                fmpz_mpoly_zero(poly)
            elif x == 1:
                fmpz_mpoly_one(poly)
            else:
                cdef fmpz_t c
                fmpz_init(c)
                if isinstance(x, Integer):
                    fmpz_set_mpz(c, x.value)
                else:
                    fmpz_set_si(c, x)
                fmpz_mpoly_set_coeff_fmpz(poly, c, NULL)
                fmpz_clear(c)
            return self._new_element(poly)

        elif isinstance(x, str):
            cdef fmpz_mpoly_t poly
            fmpz_mpoly_init(poly, self._n, self._flint_order)
            if fmpz_mpoly_set_str(poly, x.encode('utf-8')) != 0:
                raise ValueError(f"Could not parse polynomial from string: {x}")
            return self._new_element(poly)

        elif isinstance(x, MPolynomial_integer_dense_flint):
            return x

        else:
            raise TypeError(f"Cannot convert {type(x)} to a multivariate polynomial")


cdef class MPolynomial_integer_dense_flint(MPolynomialElement):
    """
    An element of a multivariate polynomial ring over the integers using FLINT's ``fmpz_mpoly``.

    EXAMPLES::

        sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
        sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
        sage: x, y = R.gens()
        sage: p = x^2 + y
        sage: p
        x^2 + y
    """
    def __init__(self, parent, poly):
        """
        Initialize a multivariate polynomial over ZZ using FLINT.

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: p = x + y
            sage: p._poly is not None  # indirect doctest
            True
        """
        self._parent = parent
        self._poly = poly

    def __dealloc__(self):
        """
        Free the memory allocated for the FLINT polynomial.
        """
        if self._poly is not None:
            fmpz_mpoly_clear(self._poly)

    def __repr__(self):
        """
        Return a string representation of this polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: repr(x + y)
            'x + y'
        """
        return fmpz_mpoly_get_str(self._poly).decode('utf-8')

    def __str__(self):
        """
        Return a string representation of this polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: str(x + y)
            'x + y'
        """
        return fmpz_mpoly_get_str(self._poly).decode('utf-8')

    def __hash__(self):
        """
        Return a hash of this polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: hash(x + y) == hash(x + y)
            True
        """
        return hash(str(self))

    def __add__(self, other):
        """
        Add two polynomials.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x^2 + y) + (x*y + 1)
            x^2 + x*y + y + 1
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_add(result, self._poly, other._poly)
        return self._parent._new_element(result)

    def __sub__(self, other):
        """
        Subtract two polynomials.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x^2 + y) - (x*y + 1)
            x^2 - x*y + y - 1
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_sub(result, self._poly, other._poly)
        return self._parent._new_element(result)

    def __mul__(self, other):
        """
        Multiply two polynomials or by a scalar.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x + y) * (x - y)
            x^2 - y^2
            sage: (x + y) * 2
            2*x + 2*y
            sage: 2 * (x + y)
            2*x + 2*y
        """
        if isinstance(other, (int, Integer)):
            cdef fmpz_t c
            cdef fmpz_mpoly_t result
            fmpz_init(c)
            if isinstance(other, Integer):
                fmpz_set_mpz(c, other.value)
            else:
                fmpz_set_si(c, other)
            fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
            fmpz_mpoly_scalar_mul_fmpz(result, self._poly, c)
            fmpz_clear(c)
            return self._parent._new_element(result)
        elif not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_mul(result, self._poly, other._poly)
        return self._parent._new_element(result)

    def __rmul__(self, other):
        """
        Right multiplication by a scalar.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: 3 * (x + y)
            3*x + 3*y
        """
        return self.__mul__(other)

    def __neg__(self):
        """
        Negate the polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: -(x + y)
            -x - y
        """
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_neg(result, self._poly)
        return self._parent._new_element(result)

    def __pow__(self, exponent):
        """
        Raise the polynomial to a power.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x + y)^2
            x^2 + 2*x*y + y^2
            sage: (x + y)^0
            1
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
        return self._parent._new_element(result)

    def __eq__(self, other):
        """
        Check equality with another polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x + y) == (y + x)
            True
            sage: (x + y) == (x + y + 1)
            False
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            return False
        return bool(fmpz_mpoly_equal(self._poly, other._poly))

    def __ne__(self, other):
        """
        Check inequality with another polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x + y) != (y + x)
            False
            sage: (x + y) != (x + y + 1)
            True
        """
        return not self.__eq__(other)

    # Properties
    def is_zero(self):
        """
        Check if the polynomial is zero.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: R.zero().is_zero()
            True
            sage: (x + y).is_zero()
            False
        """
        return bool(fmpz_mpoly_is_zero(self._poly))

    def is_one(self):
        """
        Check if the polynomial is one.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: R.one().is_one()
            True
            sage: (x + y).is_one()
            False
        """
        return bool(fmpz_mpoly_is_one(self._poly))

    def total_degree(self):
        """
        Return the total degree of the polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x^2 * y + x * y^3).total_degree()
            4
        """
        return fmpz_mpoly_total_degree(self._poly)

    def monomials(self):
        """
        Return the list of monomials of this polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x^2 * y + x).monomials()
            [(2, 1), (1, 0)]
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
        Return the list of coefficients of this polynomial as Sage Integers.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (2*x^2 * y + 3*x).coefficients()
            [3, 2]
        """
        cdef slong i, length = fmpz_mpoly_length(self._poly)
        cdef fmpz_t coeff
        cdef mpz_t mpz_coeff
        cdef slong *exps = <slong*>malloc(self._parent._n * sizeof(slong))
        coefficients = []
        fmpz_init(coeff)
        for i in range(length):
            fmpz_mpoly_get_term_coeff_fmpz(coeff, self._poly, i)
            mpz_coeff = fmpz_get_mpz(coeff)
            coefficients.append(Integer.from_mpz_t(mpz_coeff))
        fmpz_clear(coeff)
        free(exps)
        return coefficients

    def subs(self, **kwargs):
        """
        Substitute variables with values.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x^2 + y).subs(x=1, y=2)
            3
            sage: (x^2 + y).subs(x=Integer(1), y=Integer(2))
            3
        """
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_set(result, self._poly)

        for var, value in kwargs.items():
            if not isinstance(value, (int, Integer)):
                raise TypeError(f"Substitution value for {var} must be an integer or Integer")
            try:
                idx = self._parent._names.index(var)
            except ValueError:
                raise ValueError(f"Variable {var} not in the ring")

            cdef fmpz_t val
            fmpz_init(val)
            if isinstance(value, Integer):
                fmpz_set_mpz(val, value.value)
            else:
                fmpz_set_si(val, value)
            fmpz_mpoly_evaluate_one_fmpz(result, result, idx, val)
            fmpz_clear(val)

        return self._parent._new_element(result)

    def gcd(self, other):
        """
        Compute the GCD of this polynomial and another.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_integer_dense_flint import MPolynomialRing_integer_dense_flint
            sage: R = MPolynomialRing_integer_dense_flint(2, ['x', 'y'])
            sage: x, y = R.gens()
            sage: (x^2 - y^2).gcd(x - y)
            x - y
        """
        if not isinstance(other, MPolynomial_integer_dense_flint):
            other = self._parent(other)
        cdef fmpz_mpoly_t result
        fmpz_mpoly_init(result, self._parent._n, self._parent._flint_order)
        fmpz_mpoly_gcd(result, self._poly, other._poly)
        return self._parent._new_element(result)
