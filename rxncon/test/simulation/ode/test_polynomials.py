import pytest

from rxncon.simulation.ode.polynomials import Polynomial, PolynomialTerm, Symbol, TrivialMonomial, Monomial, \
    TrivialMulSymbol, MonomialFactor


### MONOMIAL ALGEBRA ###
def test_monomial_multiplication_positive_powers() -> None:
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    first_monomial = Monomial({MonomialFactor(x, 2), MonomialFactor(y, 3), MonomialFactor(z, 4)})
    second_monomial = Monomial({MonomialFactor(x, 1), MonomialFactor(y, 1), MonomialFactor(z, 1)})

    assert first_monomial * second_monomial == \
           Monomial({MonomialFactor(x, 3), MonomialFactor(y, 4), MonomialFactor(z, 5)})


def test_monomial_multiplication_mixed_sign_powers() -> None:
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    first_monomial = Monomial({MonomialFactor(x, 2), MonomialFactor(y, 3), MonomialFactor(z, 4)})
    second_monomial = Monomial({MonomialFactor(x, 1), MonomialFactor(y, -3), MonomialFactor(z, 1)})

    assert first_monomial * second_monomial == \
           Monomial({MonomialFactor(x, 3), MonomialFactor(z, 5)})


def test_monomial_multiplication_to_trivial() -> None:
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    first_monomial = Monomial({MonomialFactor(x, 2), MonomialFactor(y, -3), MonomialFactor(z, 4)})
    second_monomial = Monomial({MonomialFactor(x, -2), MonomialFactor(z, -4), MonomialFactor(y, 3)})

    assert first_monomial * second_monomial == Monomial({MonomialFactor(TrivialMulSymbol(), 0)})
    assert second_monomial * first_monomial == TrivialMonomial()


def test_monomial_symbols() -> None:
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    Monomial({MonomialFactor(x, 2), MonomialFactor(y, 3), MonomialFactor(z, 4)})

    assert Monomial({MonomialFactor(x, 2), MonomialFactor(y, 3), MonomialFactor(z, 4)}).symbols == \
           {x, y, z}

    assert Monomial({MonomialFactor(x, 2), MonomialFactor(y, 3)}).symbols == {x, y}


### POLYNOMIAL ALGEBRA ###
def test_polynomial_multiplication_by_int(x_plus_y: Polynomial, two_x_plus_two_y: Polynomial) -> None:
    assert x_plus_y * 2 == two_x_plus_two_y
    assert 2 * x_plus_y == two_x_plus_two_y


def test_polynomial_addition(x_plus_y: Polynomial, x_minus_y: Polynomial, two_x: Polynomial, two_x_plus_two_y: Polynomial,
                             two_x_minus_two_y: Polynomial) -> None:
    assert x_plus_y + x_plus_y == two_x_plus_two_y
    assert x_plus_y + x_minus_y == two_x
    assert x_minus_y + x_minus_y == two_x_minus_two_y


def test_polynomial_addition_by_scalar(one_x: Polynomial, two_x: Polynomial, const_one: Polynomial) -> None:
    assert 2 * (one_x + 1) == 2 * (one_x + const_one) == two_x + 2 == two_x + 2 * const_one == two_x + const_one * 2


def test_polynomial_subtraction(one_x: Polynomial, one_y: Polynomial, two_x_minus_two_y: Polynomial) -> None:
    assert 2 * one_x - 2 * one_y == two_x_minus_two_y


def test_polynomial_symbols(two_x: Polynomial, x_sq_minus_y_sq: Polynomial) -> None:
    x = Symbol('x')
    y = Symbol('y')

    assert two_x.symbols == {x}
    assert x_sq_minus_y_sq.symbols == {x, y}


### TEST FIXTURES ###
@pytest.fixture
def const_one() -> Polynomial:
    x = TrivialMulSymbol()
    one = PolynomialTerm(Monomial({MonomialFactor(x, 0)}), 1.0)

    return Polynomial({one})

@pytest.fixture
def one_x() -> Polynomial:
    x = Symbol('x')
    one_x = PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 1.0)

    return Polynomial({one_x})


@pytest.fixture
def one_y() -> Polynomial:
    y = Symbol('y')
    one_y = PolynomialTerm(Monomial({MonomialFactor(y, 1)}), 1.0)

    return Polynomial({one_y})


@pytest.fixture
def two_x() -> Polynomial:
    x = Symbol('x')
    term_two_x = PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 2.0)

    return Polynomial({term_two_x})


@pytest.fixture
def x_plus_y() -> Polynomial:
    x = Symbol('x')
    y = Symbol('y')
    term_x = PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 1.0)
    term_y = PolynomialTerm(Monomial({MonomialFactor(y, 1)}), 1.0)

    return Polynomial({term_x, term_y})


@pytest.fixture
def two_x_plus_two_y() -> Polynomial:
    x = Symbol('x')
    y = Symbol('y')
    term_two_x = PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 2.0)
    term_two_y = PolynomialTerm(Monomial({MonomialFactor(y, 1)}), 2.0)

    return Polynomial({term_two_x, term_two_y})


@pytest.fixture
def two_x_minus_two_y() -> Polynomial:
    x = Symbol('x')
    y = Symbol('y')
    term_two_x = PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 2.0)
    term_minus_two_y = PolynomialTerm(Monomial({MonomialFactor(y, 1)}), -2.0)

    return Polynomial({term_two_x, term_minus_two_y})


@pytest.fixture
def x_minus_y() -> Polynomial:
    x = Symbol('x')
    y = Symbol('y')
    term_x = PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 1.0)
    term_minus_y = PolynomialTerm(Monomial({MonomialFactor(y, 1)}), -1.0)

    return Polynomial({term_x, term_minus_y})


@pytest.fixture
def x_sq_minus_y_sq() -> Polynomial:
    x = Symbol('x')
    y = Symbol('y')
    term_x_sq = PolynomialTerm(Monomial({MonomialFactor(x, 2)}), 1.0)
    term_minus_y_sq = PolynomialTerm(Monomial({MonomialFactor(y, 2)}), -1.0)

    return Polynomial({term_x_sq, term_minus_y_sq})


@pytest.fixture
def x_sq_plus_two_xy_plus_y_sq() -> Polynomial:
    x = Symbol('x')
    y = Symbol('y')
    term_2xy = PolynomialTerm(Monomial({MonomialFactor(x, 1), MonomialFactor(y, 1)}), 2.0)
    term_x_sq = PolynomialTerm(Monomial({MonomialFactor(x, 2)}), 1.0)
    term_y_sq = PolynomialTerm(Monomial({MonomialFactor(y, 2)}), 1.0)

    return Polynomial({term_x_sq, term_y_sq, term_2xy})