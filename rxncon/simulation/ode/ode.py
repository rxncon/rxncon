import typecheck as tc

import rxncon.simulation.ode.polynomials as pol


class ODE:
    def __init__(self, time_derivative: pol.Symbol, polynomial: pol.Polynomial):
        self.time_derivative = time_derivative
        self.polynomial = polynomial

