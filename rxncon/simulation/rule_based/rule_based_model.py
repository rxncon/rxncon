from typing import Set, Optional
from typecheck import typecheck

from rxncon.semantics.molecule import MoleculeDefinition, Molecule
from rxncon.semantics.rule import Rule

class RuleBasedModel:
    @typecheck
    def __init__(self, molecule_defs: Set[MoleculeDefinition], rules: Set['Rule'],
                 parameters: Set['Parameter'], initial_conditions: Set['InitialCondition']):
        self.mol_defs = molecule_defs
        self.rules = rules
        self.parameters = parameters
        self.initial_conditions = initial_conditions

        self._validate()

    def set_parameter_value(self, parameter_name, parameter_value):
        # @todo
        pass

    def set_initial_condition(self, molecule, value):
        # @todo
        pass

    def _validate(self):
        for initial_condition in self.initial_conditions:
            if initial_condition.molecule.mol_def not in self.mol_defs:
                raise ValueError('Initial condition {0} refers to unknown molecule def {1}.'
                                 .format(initial_condition, initial_condition.molecule.mol_def))

        for rule in self.rules:
            for molecule in rule.molecules:
                if molecule not in self.mol_defs:
                    raise ValueError('Rule {0} contains molecule def {1}, which is absent in the model'
                                     .format(rule, molecule))


class Parameter:
    @typecheck
    def __init__(self, name: str, value: Optional[str]):
        self.name = name
        self.value = value

    @typecheck
    def __eq__(self, other: 'Parameter') -> bool:
        return self.name == other.name and self.value == other.value

    def __hash__(self):
        return hash(self.name)

    def __lt__(self, other: 'Parameter'):
        if self.name != other.name:
            return self.name < other.name
        if not self.value:
            return True
        elif not other.value:
            return False
        else:
            return self.value < other.value

    def __repr__(self):
        return str(self)

    def __str__(self) -> bool:
        return 'Parameter: {0} = {1}'.format(self.name, self.value)


class InitialCondition:
    def __init__(self, molecule: Molecule, value):
        self.molecule = molecule
        self.value = value

    def __eq__(self, other: 'InitialCondition'):
        return self.molecule == other.molecule and self.value == other.value

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'InitialCondition: {0} = {1}'.format(self.molecule, self.value)
