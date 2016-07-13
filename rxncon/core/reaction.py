from typing import Dict, Tuple, List, Optional
import re
from typecheck import typecheck

from rxncon.core.spec import MolSpec, MRnaSpec, ProteinSpec, LocusResolution, DnaSpec, mol_spec_from_string
from rxncon.core.state import StateDef, State, state_from_string, STATE_DEFS
from rxncon.util.utils import members

class Reactant:
    @typecheck
    def __init__(self, component: MolSpec, state: Optional[State]):
        assert component.is_component_spec
        self.component, self.state = component, state

    def __str__(self) -> str:
        return 'Reactant<{0}>:{1}'.format(str(self.component), str(self.state))

    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def __eq__(self, other: 'Reactant'):
        return self.component == other.component and self.state == other.state


class ReactionDef:
    ARROW_TWO_HEADS = '<->'
    ARROW_ONE_HEAD = '->'
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'

    @typecheck
    def __init__(self, state_defs: List[StateDef], name: str, representation_def: str, variables_def: Dict[str, Tuple],
                 reactants_def: str):
        self.name, self.state_defs, self.representation_def, self.variables_def, self.reactants_defs = \
            name, state_defs, representation_def, variables_def, reactants_def

        self._parse_reactants_def()

    @typecheck
    def __eq__(self, other: 'ReactionDef') -> bool:
        return self.state_defs == other.state_defs and self.name == other.name and self.representation_def == other.representation_def \
            and self.variables_def == other.variables_def and self.reactants_defs == other.reactants_defs

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'ReactionDef: {0}; representation_def: {1}; reactants_defs: {2} '\
            .format(self.name, self.representation_def, self.reactants_defs)

    def _parse_reactants_def(self):
        if self.ARROW_TWO_HEADS in self.reactants_defs:
            arrow = self.ARROW_TWO_HEADS
        elif self.ARROW_ONE_HEAD in self.reactants_defs:
            arrow = self.ARROW_ONE_HEAD
        else:
            raise AssertionError('Reaction definition requires presence of an arrow')

        reactants_def_pre_str, reactants_def_post_str = self.reactants_defs.split(arrow)

        self.reactant_defs_pre = [x.strip() for x in reactants_def_pre_str.split('+')]
        self.reactant_defs_post = [x.strip() for x in reactants_def_post_str.split('+')]

    def _parse_reactant(self, definition: str, index, variables: Dict) -> Reactant:
        def parse_state(state_str: str) -> State:
            for var_symbol, val in variables.items():
                for method in members(val):
                    var_with_method = '{0}.{1}'.format(var_symbol, method)
                    if var_with_method in state_str:
                        state_str = state_str.replace(var_with_method, str(getattr(val, method)()))

                if var_symbol in state_str:
                    state_str = state_str.replace(var_symbol, str(val))

            return state_from_string(self.state_defs, state_str)

        def parse_component(component_str: str) -> MolSpec:
            component_parts = component_str.split('.')
            assert len(component_parts) < 3

            component = variables[component_parts[0]].to_component_spec()
            method_str = None if len(component_parts) == 1 else component_parts[1]

            if method_str:
                try:
                    method = getattr(component, method_str)
                except AttributeError:
                    raise SyntaxError('Syntax error: {}'.format(component_str))

                # Some of the 'method calls' are actually properties.
                component = method if not callable(method) else method()

            return component

        component_str, state_str = definition.split('#')

        component, state = parse_component(component_str), parse_state(state_str)
        component.structure_index = index

        return Reactant(component, state)

    @typecheck
    def matches_representation(self, representation: str) -> bool:
        return True if re.match(self._to_matching_regex(), representation) else False

    @typecheck
    def representation_from_variables(self, variables: Dict) -> str:
        representation = self.representation_def
        for var, val in variables.items():
            representation = representation.replace(var, str(val))

        return representation

    @typecheck
    def variables_from_representation(self, representation: str) -> Dict:
        assert self.matches_representation(representation)

        variables = {}
        for var, var_def in self.variables_def.items():
            var_regex = self._to_base_regex().replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            val_spec = mol_spec_from_string(val_str)
            assert isinstance(val_spec, self.variables_def[var][0]), \
                '{0} is of type {1}, required to be of type {2}'.format(var, type(val_spec), self.variables_def[var][0])
            assert val_spec.has_resolution(self.variables_def[var][1]), \
                '{0} is of resolution {1}, required to be of resolution {2}'.format(var, val_spec.resolution, self.variables_def[var][1])

            variables[var] = val_spec

        return variables

    @typecheck
    def reactants_pre_from_variables(self, variables: Dict) -> List[Reactant]:
        return [self._parse_reactant(x, i, variables) for i, x in enumerate(self.reactant_defs_pre)]

    @typecheck
    def reactants_post_from_variables(self, variables: Dict) -> List[Reactant]:
        return [self._parse_reactant(x, i, variables) for i, x in enumerate(self.reactant_defs_post)]

    @typecheck
    def _to_base_regex(self) -> str:
        return '^{}$'.format(self.representation_def.replace('+', '\+'))

    @typecheck
    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.variables_def.keys():
            regex = regex.replace(var, self.SPEC_REGEX_GROUPED)

        return regex


REACTION_DEFINITIONS = [
    ReactionDef(
        STATE_DEFS,
        'phosphorylation',
        '$x_p+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x# + $y#$y-{0} -> $x# + $y#$y-{p}'
    ),
    ReactionDef(
        STATE_DEFS,
        'phosphotransfer',
        '$x_pt_$y',
        {
            '$x': (ProteinSpec, LocusResolution.residue),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}'
    ),
    ReactionDef(
        STATE_DEFS,
        'protein-protein-interaction',
        '$x_ppi_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'transcription',
        '$x_trsc_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (DnaSpec, LocusResolution.component)
        },
        '$x# + $y# -> $x# + $y# + $y.to_mrna_component_spec#'
    ),
    ReactionDef(
        STATE_DEFS,
        'translation',
        '$x_trsl_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (MRnaSpec, LocusResolution.component)
        },
        '$x# + $y# -> $x# + $y# + $y.to_protein_component_spec#'
    ),
    ReactionDef(
        STATE_DEFS,
        'intra-protein-interaction',
        '$x_ipi_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x--0,$y--0 -> $x#$x--[$y.locus]'
    ),
    ReactionDef(
        STATE_DEFS,
        'gene-protein-interaction',
        '$x_bind_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (DnaSpec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x#$x--$y + $y#$x--$y'
    )
]


class Reaction:
    @typecheck
    def __init__(self, definition: ReactionDef, variables: Dict):
        self.definition, self.variables = definition, variables

    def __str__(self) -> str:
        return self.definition.representation_from_variables(self.variables)

    @typecheck
    def __eq__(self, other: 'Reaction') -> bool:
        return self.definition == other.definition and self.variables == other.variables

    @property
    def reactants_pre(self) -> List[Reactant]:
        return self.definition.reactants_pre_from_variables(self.variables)

    @property
    def reactants_post(self) -> List[Reactant]:
        return self.definition.reactants_post_from_variables(self.variables)

    @property
    def sources(self) -> List[State]:
        return [reactant.state for reactant in self.reactants_pre]

    @property
    def products(self) -> List[State]:
        return [reactant.state for reactant in self.reactants_post]


def reaction_from_string(reaction_defs: List[ReactionDef], representation: str) -> Reaction:
    the_definition = next((reaction_def for reaction_def in reaction_defs
                          if reaction_def.matches_representation(representation)), None)

    assert the_definition, 'Could not match reaction {} with definition'.format(representation)
    variables = the_definition.variables_from_representation(representation)

    return Reaction(the_definition, variables)
