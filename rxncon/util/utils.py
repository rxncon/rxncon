import functools
from enum import Enum

def compose(*functions):
    return functools.reduce(lambda f, g: lambda x: f(g(x)), functions, lambda x: x)

class OrderedEnum(Enum):

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0}: {1}".format(self.name, self.value)

    def __lt__(self, other):
        return self.value < other.value