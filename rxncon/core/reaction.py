from typing import Dict, Any, List, Optional, Union
import re
from typecheck import typecheck

from rxncon.core.spec import Spec, BondSpec, MolSpec, MRNASpec, ProteinSpec, LocusResolution, DNASpec, mol_spec_from_string, spec_from_string
from rxncon.core.state import StateDef, State, state_from_string, STATE_DEFS, FullyNeutralState
from rxncon.util.utils import members


class ReactionTerm:
    pass


class MoleculeReactionTerm(ReactionTerm):
    @typecheck
    def __init__(self, spec: MolSpec, states: List[State], bonds: List[BondSpec]):
        assert spec.is_component_spec
        self.spec, self.states, self.bonds = spec, states, bonds

    @typecheck
    def __eq__(self, other: ReactionTerm) -> bool:
        return isinstance(other, MoleculeReactionTerm) and self.states == other.states and self.bonds == other.bonds

    def __str__(self) -> str:
        return 'MoleculeReactionTerm<{0}>states:{1},bonds:{2}'\
            .format(str(self.spec), ''.join(str(x) for x in self.states), ''.join(str(x) for x in self.bonds))

    def __repr__(self) -> str:
        return str(self)

    @property
    def is_fully_neutral(self) -> bool:
        return FullyNeutralState() in self.states and len(self.states) == 1


class BondReactionTerm(ReactionTerm):
    @typecheck
    def __init__(self, spec: BondSpec, states: List[State]):
        self.spec, self.states = spec, states

    @typecheck
    def __eq__(self, other: ReactionTerm) -> bool:
        return isinstance(other, BondReactionTerm) and self.spec == other.spec and self.states == other.states

    def __str__(self) -> str:
        return 'BondReactionTerm<{0}>states:{1}'.format(str(self.spec), ''.join(str(x) for x in self.states))

    def __repr__(self) -> str:
        return str(self)


class ReactionDef:
    ARROW_TWO_HEADS = '<->'
    ARROW_ONE_HEAD = '->'
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'

    @typecheck
    def __init__(self, state_defs: List[StateDef], name: str, representation_def: str, variables_def: Dict[str, Any],
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

    @typecheck
    def matches_representation(self, representation: str) -> bool:
        return True if re.match(self._to_matching_regex(), representation) else False

    @typecheck
    def representation_from_variables(self, variables: Dict[str, Any]) -> str:
        representation = self.representation_def
        for var, val in variables.items():
            representation = representation.replace(var, str(val))

        return representation

    @typecheck
    def variables_from_representation(self, representation: str) -> Dict[str, Any]:
        assert self.matches_representation(representation)

        variables = {}
        for var, var_def in self.variables_def.items():
            var_regex = self._to_base_regex().replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            val_spec = mol_spec_from_string(val_str)

            variables[var] = val_spec

        return variables

    @typecheck
    def validate_variables(self, variables: Dict[str, Any]):
        for var, val in variables.items():
            assert isinstance(val, self.variables_def[var][0]), \
                '{0} is of type {1}, required to be of type {2}'.format(var, type(val), self.variables_def[var][0])
            assert val.has_resolution(self.variables_def[var][1]), \
                '{0} is of resolution {1}, required to be of resolution {2}'.format(var, val.resolution, self.variables_def[var][1])

    @typecheck
    def terms_lhs_from_variables(self, variables: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, i, variables) for i, x in enumerate(self.reactant_defs_lhs)]

    @typecheck
    def terms_rhs_from_variables(self, variables: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, i, variables) for i, x in enumerate(self.reactant_defs_rhs)]

    def _parse_reactants_def(self):
        if self.ARROW_TWO_HEADS in self.reactants_defs:
            arrow = self.ARROW_TWO_HEADS
        elif self.ARROW_ONE_HEAD in self.reactants_defs:
            arrow = self.ARROW_ONE_HEAD
        else:
            raise AssertionError('Reaction definition requires presence of an arrow')

        reactants_def_lhs_str, reactants_def_rhs_str = self.reactants_defs.split(arrow)

        self.reactant_defs_lhs = [x.strip() for x in reactants_def_lhs_str.split('+')]
        self.reactant_defs_rhs = [x.strip() for x in reactants_def_rhs_str.split('+')]

    def _parse_term(self, definition: str, index, variables: Dict[str, Any]) -> ReactionTerm:
        def parse_bond_spec_values(values_str: str):
            if not values_str:
                return []

            try:
                potential_spec_str = values_str
                for var_symb, var_val in variables.items():
                    potential_spec_str = potential_spec_str.replace(var_symb, str(var_val))

                spec = spec_from_string(potential_spec_str)
                if isinstance(spec, BondSpec):
                    return [spec]
                else:
                    return []
            except SyntaxError:
                return []

        def parse_state_values(values_str: str):
            if not values_str:
                return []

            try:
                return [parse_state(value) for value in values_str.split(',')]
            except SyntaxError:
                return []

        def parse_state(state_str: str) -> State:
            for var_symbol, var_val in variables.items():
                for method in members(var_val):
                    var_with_method = '{0}.{1}'.format(var_symbol, method)
                    if var_with_method in state_str:
                        method_res = getattr(var_val, method)
                        method_res = method_res if not callable(method_res) else method_res()
                        state_str = state_str.replace(var_with_method, str(method_res))

                if var_symbol in state_str:
                    state_str = state_str.replace(var_symbol, str(var_val))

            return state_from_string(state_str)

        def parse_molecule_or_bond(spec_str: str) -> Spec:
            if '~' in spec_str:
                for var_symb, var_val in variables.items():
                    spec_str = spec_str.replace(var_symb, str(var_val))

                return spec_from_string(spec_str)

            component_parts = spec_str.split('.')
            assert len(component_parts) < 3

            component = variables[component_parts[0]].to_component_spec()
            method_str = None if len(component_parts) == 1 else component_parts[1]

            if method_str:
                try:
                    method = getattr(component, method_str)
                except AttributeError:
                    raise SyntaxError('Syntax error: {}'.format(spec_str))

                # Some of the 'method calls' are actually properties.
                component = method if not callable(method) else method()

            return component

        spec_str, values_str = definition.split('#')

        spec = parse_molecule_or_bond(spec_str)
        spec.structure_index = index

        if isinstance(spec, BondSpec):
            return BondReactionTerm(spec, parse_state_values(values_str))
        elif isinstance(spec, MolSpec):
            return MoleculeReactionTerm(spec, parse_state_values(values_str), parse_bond_spec_values(values_str))


    @typecheck
    def _to_base_regex(self) -> str:
        return '^{}$'.format(self.representation_def.replace('+', '\+'))

    @typecheck
    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.variables_def.keys():
            regex = regex.replace(var, self.SPEC_REGEX_GROUPED)

        return regex


REACTION_DEFS = [
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
        'dephosphorylation',
        '$x_p-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x# + $y#$y-{p} -> $x# + $y#$y-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'auto-phosphorylation',
        '$x_ap+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{0} -> $y#$y-{p}'
    ),
    ReactionDef(
        STATE_DEFS,
        'ubiquitination',
        '$x_ub+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x# + $y#$y-{0} -> $x# + $y#$y-{ub}'
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
        '$x_ppi+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x#$x~$y + $y#$x~$y + $x~$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'protein-protein-dissociation',
        '$x_ppi-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x~$y + $y#$x~$y + $x~$y#$x--$y -> $x#$x--0 + $y#$y--0'
    ),
    ReactionDef(
        STATE_DEFS,
        'interaction',
        '$x_i+_$y',
        {
            '$x': (MolSpec, LocusResolution.domain),
            '$y': (MolSpec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x#$x~$y + $y#$x~$y + $x~$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'dissociation',
        '$x_i-_$y',
        {
            '$x': (MolSpec, LocusResolution.domain),
            '$y': (MolSpec, LocusResolution.domain)
        },
        '$x#$x~$y + $y#$x~$y + $x~$y#$x--$y -> $x#$x--0 + $y#$y--0'
    ),
    ReactionDef(
        STATE_DEFS,
        'transcription',
        '$x_trsc_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (DNASpec, LocusResolution.component)
        },
        '$x# + $y# -> $x# + $y# + $y.to_mrna_component_spec#0'
    ),
    ReactionDef(
        STATE_DEFS,
        'translation',
        '$x_trsl_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (MRNASpec, LocusResolution.component)
        },
        '$x# + $y# -> $x# + $y# + $y.to_protein_component_spec#0'
    ),
    ReactionDef(
        STATE_DEFS,
        'intra-protein-interaction',
        '$x_ipi_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x--0,$y--0 <-> $x#$x~$y + $x~$y#$x--[$y.locus]'
    ),
    ReactionDef(
        STATE_DEFS,
        'gene-protein-interaction',
        '$x_bind_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (DNASpec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x#$x~$y + $y#$x~$y + $x~$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'protein-degradation',
        '$x_deg_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$x# + $y# -> $x#'
    )
]


class Reaction:
    @typecheck
    def __init__(self, definition: ReactionDef, variables: Dict[str, Any]):
        self.name            = definition.name
        self.terms_lhs       = definition.terms_lhs_from_variables(variables)
        self.terms_rhs       = definition.terms_rhs_from_variables(variables)
        self._representation = definition.representation_from_variables(variables)

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self._representation

    @typecheck
    def __eq__(self, other: 'Reaction') -> bool:
        return self.terms_lhs == other.terms_lhs and self.terms_rhs == other.terms_rhs and str(self) == str(other)

    @property
    @typecheck
    def components_lhs(self) -> List[MolSpec]:
        return [x.spec for x in self.terms_lhs if isinstance(x, MoleculeReactionTerm)]

    @property
    @typecheck
    def components_rhs(self) -> List[MolSpec]:
        return [x.spec for x in self.terms_rhs if isinstance(x, MoleculeReactionTerm)]

    @property
    @typecheck
    def bonds_lhs(self) -> List[BondSpec]:
        return [x.spec for x in self.terms_lhs if isinstance(x, BondReactionTerm)]

    @property
    @typecheck
    def bonds_rhs(self) -> List[BondSpec]:
        return [x.spec for x in self.terms_rhs if isinstance(x, BondReactionTerm)]

    @property
    @typecheck
    def consumed_states(self) -> List[State]:
        states = []

        for term in self.terms_lhs:
            if isinstance(term, MoleculeReactionTerm) and term.spec in self.components_rhs:
                states += term.states
            elif isinstance(term, BondReactionTerm) and term.spec.first.to_component_spec() in self.components_rhs and \
                    term.spec.second.to_component_spec() in self.components_rhs:
                states += term.states

        return states

    @property
    @typecheck
    def produced_states(self) -> List[State]:
        states = []

        for term in self.terms_rhs:
            if isinstance(term, MoleculeReactionTerm) and term.spec in self.components_lhs:
                states += term.states
            elif isinstance(term, BondReactionTerm) and term.spec.first.to_component_spec() in self.components_lhs and \
                    term.spec.second.to_component_spec() in self.components_lhs:
                states += term.states

        return states

    @property
    @typecheck
    def degraded_states(self) -> List[State]:
        states = []

        for term in self.terms_lhs:
            if isinstance(term, MoleculeReactionTerm) and term.spec not in self.components_rhs:
                states += term.states
            elif isinstance(term, BondReactionTerm) and term.spec.first.to_component_spec() not in self.components_rhs and \
                    term.spec.second.to_component_spec() not in self.components_rhs:
                states += term.states

        return states

    @property
    @typecheck
    def synthesised_states(self) -> List[State]:
        states = []

        for term in self.terms_rhs:
            if isinstance(term, MoleculeReactionTerm) and term.spec not in self.components_lhs:
                states += term.states
            elif isinstance(term, BondReactionTerm) and term.spec.first.to_component_spec() not in self.components_lhs and \
                    term.spec.second.to_component_spec() not in self.components_lhs:
                states += term.states

        return states

    @property
    @typecheck
    def degraded_components(self) -> List[MolSpec]:
        return [component for component in self.components_lhs if component not in self.components_rhs]

    @property
    @typecheck
    def synthesised_components(self) -> List[MolSpec]:
        return [component for component in self.components_rhs if component not in self.components_lhs]


@typecheck
def matching_reaction_def(representation: str) -> Optional[ReactionDef]:
    return next((reaction_def for reaction_def in REACTION_DEFS if reaction_def.matches_representation(representation)), None)


@typecheck
def reaction_from_string(representation: str, standardize=True) -> Reaction:
    def fixed_spec_types(rxn_def: ReactionDef, variables: Dict[str, Any]) -> Dict[str, Any]:
        keys = variables.keys()
        assert len(list(keys)) == 2

        for key in keys:
            required_type = rxn_def.variables_def[key][0]
            if not isinstance(variables[key], required_type):
                if required_type is DNASpec:
                    variables[key] = variables[key].to_dna_component_spec()
                elif required_type is MRNASpec:
                    variables[key] = variables[key].to_mrna_component_spec()
                elif required_type is ProteinSpec:
                    variables[key] = variables[key].to_protein_component_spec()
                else:
                    raise NotImplementedError

        return variables

    def fixed_resolutions(rxn_def: ReactionDef, variables: Dict[str, Any]) -> Dict[str, Any]:
        keys = variables.keys()
        assert len(list(keys)) == 2

        for key in keys:
            if rxn_def.variables_def[key][1] < variables[key].resolution:
                raise SyntaxError('Specified resolution for variable {0} higher than required {1}'
                                  .format(str(variables[key]), rxn_def.variables_def[key][1]))

            other = [x for x in keys if x != key][0]
            if not variables[key].has_resolution(rxn_def.variables_def[key][1]):
                if rxn_def.variables_def[key][1] == LocusResolution.domain:
                    variables[key].locus.domain = variables[other].component_name
                elif rxn_def.variables_def[key][1] == LocusResolution.residue:
                    variables[key].locus.residue = variables[other].component_name
                else:
                    raise NotImplementedError

        return variables

    the_definition = matching_reaction_def(representation)
    assert the_definition, 'Could not match reaction {} with definition'.format(representation)

    variables = the_definition.variables_from_representation(representation)

    if standardize:
        variables = fixed_spec_types(the_definition, variables)
        variables = fixed_resolutions(the_definition, variables)

    the_definition.validate_variables(variables)
    return Reaction(the_definition, variables)