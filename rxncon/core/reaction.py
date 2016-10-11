from typing import Dict, Any, List, Optional
import re

from rxncon.core.spec import Spec, BondSpec, Spec, MRNASpec, ProteinSpec, LocusResolution, DNASpec, spec_from_str, spec_from_str
from rxncon.core.state import StateDef, State, state_from_str, STATE_DEFS, FullyNeutralState
from rxncon.util.utils import members


class ReactionTerm:
    def __init__(self, specs: List[Spec], states: List[State]):
        assert all(spec.is_component_spec for spec in specs)
        self.specs, self.states = specs, states

    def __eq__(self, other: 'ReactionTerm') -> bool:
        return isinstance(other, ReactionTerm) and self.states == other.states

    def __str__(self) -> str:
        return 'ReactionTerm<{0}>states:{1}'.format(''.join(str(spec) for spec in self.specs), ''.join(str(x) for x in self.states))

    def __repr__(self) -> str:
        return str(self)

    @property
    def is_fully_neutral(self) -> bool:
        return self.states == [FullyNeutralState()]


class ReactionDef:
    ARROW_TWO_HEADS = '<->'
    ARROW_ONE_HEAD = '->'
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'

    def __init__(self, state_defs: List[StateDef], name: str, repr_def: str, vars_def: Dict[str, Any], rule_def: str):
        self.name, self.state_defs, self.repr_def, self.vars_def, self.rule_defs = name, state_defs, repr_def, vars_def, rule_def
        self._parse_reactants_def()

    def __eq__(self, other: 'ReactionDef') -> bool:
        return self.state_defs == other.state_defs and self.name == other.name and self.repr_def == other.repr_def \
            and self.vars_def == other.vars_def and self.rule_defs == other.rule_defs

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'ReactionDef: {0}; representation_def: {1}; reactants_defs: {2} '\
            .format(self.name, self.repr_def, self.rule_defs)

    def matches_repr(self, representation: str) -> bool:
        return True if re.match(self._to_matching_regex(), representation) else False

    def repr_from_vars(self, variables: Dict[str, Any]) -> str:
        representation = self.repr_def
        for var, val in variables.items():
            representation = representation.replace(var, str(val))

        return representation

    def vars_from_repr(self, representation: str) -> Dict[str, Any]:
        assert self.matches_repr(representation)

        variables = {}
        for var, var_def in self.vars_def.items():
            var_regex = self._to_base_regex().replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.vars_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            val_spec = spec_from_str(val_str)

            variables[var] = val_spec

        return variables

    def validate_vars(self, variables: Dict[str, Any]):
        for var, val in variables.items():
            assert isinstance(val, self.vars_def[var][0]), \
                '{0} is of type {1}, required to be of type {2}'.format(var, type(val), self.vars_def[var][0])
            assert val.has_resolution(self.vars_def[var][1]), \
                '{0} is of resolution {1}, required to be of resolution {2}'.format(var, val.resolution, self.vars_def[var][1])

    def terms_lhs_from_vars(self, variables: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, variables) for x in self.reactant_defs_lhs]

    def terms_rhs_from_vars(self, variables: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, variables) for x in self.reactant_defs_rhs]

    def _parse_reactants_def(self):
        if self.ARROW_TWO_HEADS in self.rule_defs:
            arrow = self.ARROW_TWO_HEADS
        elif self.ARROW_ONE_HEAD in self.rule_defs:
            arrow = self.ARROW_ONE_HEAD
        else:
            raise AssertionError('Reaction definition requires presence of an arrow')

        reactants_def_lhs_str, reactants_def_rhs_str = self.rule_defs.split(arrow)

        self.reactant_defs_lhs = [x.strip() for x in reactants_def_lhs_str.split('+')]
        self.reactant_defs_rhs = [x.strip() for x in reactants_def_rhs_str.split('+')]

    def _parse_term(self, definition: str, variables: Dict[str, Any]) -> ReactionTerm:
        def parse_states_str(states_str: str) -> List[State]:
            if not states_str:
                return []
            try:
                return [parse_state_str(x) for x in states_str.split(',')]
            except SyntaxError:
                return []

        def parse_state_str(state_str: str) -> State:
            for var_symbol, var_val in variables.items():
                for method in members(var_val):
                    var_with_method = '{0}.{1}'.format(var_symbol, method)
                    if var_with_method in state_str:
                        method_res = getattr(var_val, method)
                        method_res = method_res if not callable(method_res) else method_res()
                        state_str = state_str.replace(var_with_method, str(method_res))

                if var_symbol in state_str:
                    state_str = state_str.replace(var_symbol, str(var_val))

            return state_from_str(state_str)

        def parse_specs_str(specs_str: str) -> List[Spec]:
            return [spec_from_str(x) for x in specs_str.split(',')]

        specs_str, states_str = definition.split('#')

        return ReactionTerm(parse_specs_str(specs_str), parse_states_str(states_str))

    def _to_base_regex(self) -> str:
        return '^{}$'.format(self.repr_def.replace('+', '\+'))

    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.vars_def.keys():
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
            '$x': (Spec, LocusResolution.domain),
            '$y': (Spec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x#$x~$y + $y#$x~$y + $x~$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'dissociation',
        '$x_i-_$y',
        {
            '$x': (Spec, LocusResolution.domain),
            '$y': (Spec, LocusResolution.domain)
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
        'degradation',
        '$x_deg_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x# + $y# -> $x#'
    ),
    ReactionDef(
        STATE_DEFS,
        'synthesis',
        '$x_syn_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x# -> $x# + $y#0'
    )
]


class Reaction:
    def __init__(self, definition: ReactionDef, vars: Dict[str, Any]):
        self.name            = definition.name
        self.terms_lhs       = definition.terms_lhs_from_vars(vars)
        self.terms_rhs       = definition.terms_rhs_from_vars(vars)
        self._representation = definition.repr_from_vars(vars)

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self._representation

    def __eq__(self, other: 'Reaction') -> bool:
        return self.terms_lhs == other.terms_lhs and self.terms_rhs == other.terms_rhs and str(self) == str(other)

    @property
    def components_lhs(self) -> List[Spec]:
        return [spec for term in self.terms_lhs for spec in term]

    @property
    def components_rhs(self) -> List[Spec]:
        return [spec for term in self.terms_rhs for spec in term]

    @property
    def consumed_states(self) -> List[State]:
        states = []

        for term in self.terms_lhs:
            if all(spec in self.components_rhs for spec in term.specs):
                states += term.states

        return states

    @property
    def produced_states(self) -> List[State]:
        states = []

        for term in self.terms_rhs:
            if all(spec in self.components_lhs for spec in term.specs):
                states += term.states

        return states

    @property
    def degraded_states(self) -> List[State]:
        states = []

        for term in self.terms_lhs:
            if not any(spec in self.components_rhs for spec in term.specs):
                states += term.states

        return states

    @property
    def synthesised_states(self) -> List[State]:
        states = []

        for term in self.terms_rhs:
            if not any(spec in self.components_lhs for spec in term.specs):
                states += term.states

        return states

    @property
    def degraded_components(self) -> List[Spec]:
        return [component for component in self.components_lhs if component not in self.components_rhs]

    @property
    def synthesised_components(self) -> List[Spec]:
        return [component for component in self.components_rhs if component not in self.components_lhs]


def matching_reaction_def(representation: str) -> Optional[ReactionDef]:
    return next((reaction_def for reaction_def in REACTION_DEFS if reaction_def.matches_representation(representation)), None)


def reaction_from_str(representation: str, standardize=True) -> Reaction:
    def fixed_spec_types(rxn_def: ReactionDef, variables: Dict[str, Any]) -> Dict[str, Any]:
        keys = variables.keys()
        assert len(list(keys)) == 2

        for key in keys:
            required_type = rxn_def.vars_def[key][0]
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
            if rxn_def.vars_def[key][1] < variables[key].resolution:
                raise SyntaxError('Specified resolution for variable {0} higher than required {1}'
                                  .format(str(variables[key]), rxn_def.vars_def[key][1]))

            other = [x for x in keys if x != key][0]
            if not variables[key].has_resolution(rxn_def.vars_def[key][1]):
                if rxn_def.vars_def[key][1] == LocusResolution.domain:
                    variables[key].locus.domain = variables[other].component_name
                elif rxn_def.vars_def[key][1] == LocusResolution.residue:
                    variables[key].locus.residue = variables[other].component_name
                else:
                    raise NotImplementedError

        return variables

    the_definition = matching_reaction_def(representation)

    if not the_definition:
        raise SyntaxError('Could not match reaction {} with definition'.format(representation))

    variables = the_definition.vars_from_repr(representation)

    if standardize:
        variables = fixed_spec_types(the_definition, variables)
        variables = fixed_resolutions(the_definition, variables)

    the_definition.validate_vars(variables)
    return Reaction(the_definition, variables)
