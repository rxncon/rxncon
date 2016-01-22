from typing import List, Callable


class ODESystem:
    def __init__(self, odes: List['ODE']):
        self.odes = odes


class ODE:
    def __init__(self, time_varying_quantity: Quantity, terms: List['ODETerm']):
        self.time_varying_quantity = time_varying_quantity
        self.terms = terms


class Quantity:
    def __init__(self, value):
        # @todo assert value has certain properties?
        self.value = value


class ODETerm:
    def __init__(self, coefficient: float, factors: List[Quantity]):
        self.coefficient = coefficient
        self.factors = factors


class SolverInformation:
    def __init__(self, rhs_function: List[Callable], initial_conditions: List[float]):
        self.rhs_function = rhs_function
        self.initial_conditions = initial_conditions





