from enum import Enum
from typing import Set

from typecheck import typecheck

from rxncon.semantics.molecule import Molecule, Bond


class Rule:
    def __init__(self, lhs: Set['Complex'], rhs: Set['Complex'], arrow_type: 'Arrow', rates: Set['Parameter']):
        self.lhs = lhs
        self.rhs = rhs
        self.arrow_type = arrow_type
        self.rates = rates
        self._validate()

    def __eq__(self, other: 'Rule'):
        return self.lhs == other.lhs and self.rhs == other.rhs and self.arrow_type == other.arrow_type and self.rates == other.rates

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'Rule: {0} {1} {2}, {3}'.format('+'.join(str(x) for x in sorted(self.lhs)),
                                               self.arrow_type,
                                               '+'.join(str(x) for x in sorted(self.rhs)),
                                               ', '.join(str(x) for x in sorted(self.rates)))

    def _validate(self):
        if self.arrow_type == Arrow.irreversible and len(self.rates) != 1:
            raise ValueError('Rule {0} is irreversible and thus requires exactly one rate constant, {1} given'
                             .format(str(self), len(self.rates)))

        if self.arrow_type == Arrow.reversible and len(self.rates) != 2:
            raise ValueError('Rule {0} is reversible and thus requires exactly two rate constants, {1} given'
                             .format(str(self), len(self.rates)))


class Complex:
    @typecheck
    def __init__(self, molecules: Set[Molecule]):
        self.molecules = []

        for molecule in molecules:
            self.add_molecule(molecule)

        self.open_bonds = []
        self.closed_bonds = []

    @typecheck
    def __eq__(self, other: 'Complex'):
        return self.molecules == other.molecules

    def __hash__(self):
        return hash(str(self))

    def __lt__(self, other):
        return sorted(self.molecules) < sorted(other.molecules)

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'Comp<{0}>'.format(', '.join(sorted([str(x) for x in self.molecules])))

    def add_molecule(self, molecule: Molecule):
        for bond in molecule.bonds:
            if bond in self.open_bonds:
                self.open_bonds.remove(bond)
                self.closed_bonds.append(bond)
            else:
                self.open_bonds.append(bond)

        self.molecules.append(molecule)

    @property
    def is_graph_like(self):
        return len(self.open_bonds) == 0

    def make_graph_like(self):
        pass


class Arrow(Enum):
    irreversible = '->'
    reversible   = '<->'
