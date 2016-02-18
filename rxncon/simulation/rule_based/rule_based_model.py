from enum import Enum
from typing import List

import typecheck as tc

from rxncon.semantics.molecule_instance import MoleculeInstance


class RuleBasedModel:
    @tc.typecheck
    def __init__(self, molecule_defs: List['MoleculeDefinition'], rules: List['Rule'],
                 parameters: List['Parameter'], initial_conditions: List['InitialCondition']):
        self.molecule_defs = molecule_defs
        self.rules = rules
        self.parameters = parameters
        self.initial_conditions = initial_conditions

        self._validate()

    def _validate(self):
        for initial_condition in self.initial_conditions:
            if initial_condition.molecule_specification.molecule_def not in self.molecule_defs:
                raise ValueError('Initial condition {0} refers to unknown molecule def {1}.'
                                 .format(initial_condition, initial_condition.molecule_specification.molecule_def))

        for rule in self.rules:
            for molecule in rule.molecules:
                if molecule not in self.molecule_defs:
                    raise ValueError('Rule {0} contains molecule def {1}, which is absent in the model'
                                     .format(rule, molecule))


class Rule:
    @tc.typecheck
    def __init__(self, left_hand_side: List['Reactant'], right_hand_side: List['Reactant'], arrow_type: 'Arrow', rates: List['Parameter']):
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.arrow_type = arrow_type
        self.rates = rates
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Rule'):
        return self.left_hand_side == other.left_hand_side and self.right_hand_side == other.right_hand_side and \
            self.arrow_type == other.arrow_type and self.rates == other.rates

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'Rule: {0} {1} {2}, {3}'.format('+'.join(str(x) for x in self.left_hand_side), self.arrow_type,
                                               '+'.join(str(x) for x in self.right_hand_side), ', '.join(str(x) for x in self.rates))

    @property
    def molecules(self):
        molecules = []
        for side in [self.left_hand_side, self.right_hand_side]:
            if isinstance(side, MoleculeReactant) and side.molecule_specification.molecule_definition not in molecules:
                molecules.append(side.molecule_specification.molecule_def)

            elif isinstance(side, ComplexReactant):
                [molecules.append(x) for x in side.molecules if x not in molecules]

        return molecules

    def _validate(self):
        if self.arrow_type == Arrow.irreversible and len(self.rates) != 1:
            raise ValueError('Rule {0} is irreversible and thus requires exactly one rate constant, {1} given'
                             .format(str(self), len(self.rates)))

        if self.arrow_type == Arrow.reversible and len(self.rates) != 2:
            raise ValueError('Rule {0} is reversible and thus requires exactly two rate constants, {1} given'
                             .format(str(self), len(self.rates)))


class Reactant:
    pass


class MoleculeReactant(Reactant):
    @tc.typecheck
    def __init__(self, molecule_specification: MoleculeInstance):
        self.molecule_specification = molecule_specification

    @tc.typecheck
    def __eq__(self, other: Reactant) -> bool:
        return isinstance(other, MoleculeReactant) and self.molecule_specification == other.molecule_specification

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'MoleculeReactant: [{0}]'.format(self.molecule_specification)


class ComplexReactant(Reactant):
    @tc.typecheck
    def __init__(self, molecules: List[MoleculeInstance], bindings: List['Binding']):
        self.molecules = molecules
        self.bindings = bindings
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Reactant'):
        return isinstance(other, ComplexReactant) and self.molecules == other.molecules and self.bindings == other.bindings

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'ComplexReactant: Molecules = [{0}], Bindings = [{1}]'\
            .format(', '.join(str(x) for x in self.molecules), ', '.join(str(x) for x in self.bindings))

    def _validate(self):
        unique_localizations = {molecule.localization_instance for molecule in self.molecules}
        if len(unique_localizations) > 1:
            raise ValueError('Molecules making up a ComplexReactant cannot be in different localizations: {0}.'
                             .format(', '.join(str(x) for x in unique_localizations)))


class Arrow(Enum):
    irreversible = '->'
    reversible   = '<->'


class Parameter:
    @tc.typecheck
    def __init__(self, name: str, value: str):
        self.name = name
        self.value = value

    @tc.typecheck
    def __eq__(self, other: 'Parameter') -> bool:
        assert isinstance(other, Parameter)
        return self.name == other.name and self.value == other.value

    def __repr__(self):
        return str(self)

    def __str__(self) -> bool:
        return 'Parameter: {0} = {1}'.format(self.name, self.value)


class InitialCondition:
    def __init__(self, molecule_specification: MoleculeInstance, value):
        self.molecule_specification = molecule_specification
        self.value = value

    def __eq__(self, other):
        assert isinstance(other, InitialCondition)
        return self.molecule_specification == other.molecule_specification and self.value == other.value

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'InitialCondition: {0} = {1}'.format(self.molecule_specification, self.value)
