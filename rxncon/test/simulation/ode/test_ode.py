import rxncon.simulation.ode.ode as ode


### MONOMIAL ALGEBRA ###
def test_monomial_multiplication_positive_powers():
    x = ode.Symbol('x')
    y = ode.Symbol('y')
    z = ode.Symbol('z')

    first_monomial = ode.Monomial({ode.MonomialFactor(x, 2), ode.MonomialFactor(y, 3), ode.MonomialFactor(z, 4)})
    second_monomial = ode.Monomial({ode.MonomialFactor(x, 1), ode.MonomialFactor(y, 1), ode.MonomialFactor(z, 1)})

    assert first_monomial * second_monomial == \
        ode.Monomial({ode.MonomialFactor(x, 3), ode.MonomialFactor(y, 4), ode.MonomialFactor(z, 5)})


def test_monomial_multiplication_mixed_sign_powers():
    x = ode.Symbol('x')
    y = ode.Symbol('y')
    z = ode.Symbol('z')

    first_monomial = ode.Monomial({ode.MonomialFactor(x, 2), ode.MonomialFactor(y, 3), ode.MonomialFactor(z, 4)})
    second_monomial = ode.Monomial({ode.MonomialFactor(x, 1), ode.MonomialFactor(y, -3), ode.MonomialFactor(z, 1)})

    assert first_monomial * second_monomial == \
        ode.Monomial({ode.MonomialFactor(x, 3), ode.MonomialFactor(z, 5)})


def test_monomial_multiplication_to_trivial():
    x = ode.Symbol('x')
    y = ode.Symbol('y')
    z = ode.Symbol('z')

    first_monomial = ode.Monomial({ode.MonomialFactor(x, 2), ode.MonomialFactor(y, -3), ode.MonomialFactor(z, 4)})
    second_monomial = ode.Monomial({ode.MonomialFactor(x, -2), ode.MonomialFactor(z, -4), ode.MonomialFactor(y, 3)})

    assert first_monomial * second_monomial == ode.Monomial({ode.MonomialFactor(ode.TrivialSymbol(), 0)})
