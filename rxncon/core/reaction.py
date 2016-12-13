from typing import Dict, Any, List, Optional
import re

from rxncon.core.spec import Spec, MRNASpec, ProteinSpec, LocusResolution, GeneSpec, spec_from_str
from rxncon.core.state import StateDef, State, state_from_str, STATE_DEFS, FullyNeutralState
from rxncon.util.utils import members

SPEC_REGEX_MATCHING     = '([A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
SPEC_REGEX_NON_MATCHING = '(?:[A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'

BIDIRECTIONAL_VERBS = [
    'ppi', 'ipi', 'i', 'bind'
]

class ReactionTerm:
    def __init__(self, specs: List[Spec], states: List[State]):
        assert all(spec.is_component_spec for spec in specs)
        self.specs, self.states = specs, states

    def __eq__(self, other: 'ReactionTerm') -> bool:
        return isinstance(other, ReactionTerm) and self.states == other.states

    def __str__(self) -> str:
        return 'ReactionTerm<{0}>states:{1}'.format(','.join(str(spec) for spec in self.specs), ','.join(str(x) for x in self.states))

    def __repr__(self) -> str:
        return str(self)

    @property
    def is_fully_neutral(self) -> bool:
        return self.states == [FullyNeutralState()]


class ReactionDef:
    ARROW = '->'

    def __init__(self, state_defs: List[StateDef], name: str, repr_def: str, vars_def: Dict[str, Any], rule_def: str):
        self.name, self.state_defs, self.repr_def, self.vars_def, self.rule_def = name, state_defs, repr_def, vars_def, rule_def
        self._parse_reactants_def()

    def __eq__(self, other: 'ReactionDef') -> bool:
        return self.state_defs == other.state_defs and self.name == other.name and self.repr_def == other.repr_def \
            and self.vars_def == other.vars_def and self.rule_def == other.rule_def

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'ReactionDef: {0}; repr_def: {1}; rule_def: {2} '.format(self.name, self.repr_def, self.rule_def)

    def matches_repr(self, repr: str) -> bool:
        return True if re.match(self._to_matching_regex(), repr) else False

    def repr_from_vars(self, vars: Dict[str, Any]) -> str:
        representation = self.repr_def
        for var, val in vars.items():
            representation = representation.replace(var, str(val))

        return representation

    def vars_from_repr(self, repr: str) -> Dict[str, Any]:
        assert self.matches_repr(repr)

        variables = {}
        for var, var_def in self.vars_def.items():
            var_regex = self._to_base_regex().replace(var, SPEC_REGEX_MATCHING)
            for other_var in self.vars_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, SPEC_REGEX_NON_MATCHING)

            val_str = re.match(var_regex, repr).group(1)
            val_spec = spec_from_str(val_str)

            variables[var] = val_spec

        return variables

    def validate_vars(self, vars: Dict[str, Any]):
        for var, val in vars.items():
            assert isinstance(val, self.vars_def[var][0]), \
                '{0} is of type {1}, required to be of type {2}'.format(var, type(val), self.vars_def[var][0])
            assert val.has_resolution(self.vars_def[var][1]), \
                '{0} is of resolution {1}, required to be of resolution {2}'.format(var, val.resolution, self.vars_def[var][1])

    def terms_lhs_from_vars(self, vars: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, vars) for x in self.reactant_defs_lhs]

    def terms_rhs_from_vars(self, vars: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, vars) for x in self.reactant_defs_rhs]

    def _parse_reactants_def(self):
        assert self.ARROW in self.rule_def

        reactants_def_lhs_str, reactants_def_rhs_str = self.rule_def.split(self.ARROW)

        self.reactant_defs_lhs = [x.strip() for x in reactants_def_lhs_str.split('+')]
        self.reactant_defs_rhs = [x.strip() for x in reactants_def_rhs_str.split('+')]

    def _parse_term(self, definition: str, vars: Dict[str, Any]) -> ReactionTerm:
        def parse_states_str(states_str: str) -> List[State]:
            if not states_str:
                return []
            try:
                return [parse_vars_str(x, state_from_str) for x in states_str.split(',')]
            except SyntaxError:
                raise SyntaxError('Could not parse states string {}'.format(states_str))

        def parse_specs_str(specs_str: str) -> List[Spec]:
            return [parse_vars_str(x, spec_from_str).to_component_spec() for x in specs_str.split(',')]

        def parse_vars_str(vars_str: str, post_func):
            for var_symbol, var_val in vars.items():
                for method in members(var_val):
                    var_with_method = '{0}.{1}'.format(var_symbol, method)
                    if var_with_method in vars_str:
                        method_res = getattr(var_val, method)
                        method_res = method_res if not callable(method_res) else method_res()
                        vars_str = vars_str.replace(var_with_method, str(method_res))

                if var_symbol in vars_str:
                    vars_str = vars_str.replace(var_symbol, str(var_val))

            return post_func(vars_str)

        if len(definition.split('#')) != 2:
            raise SyntaxError('ReactionDef syntax error: each term should be of the form Specs#States,\n'
                              'where Specs is one or more Specs separated by commas and States is zero or\n'
                              'more States separated by commas.\nThe given definition was: {}'.format(definition))

        specs_str, states_str = definition.split('#')

        return ReactionTerm(parse_specs_str(specs_str), parse_states_str(states_str))

    def _to_base_regex(self) -> str:
        # The (?i) makes the regex case insensitive.
        return '(?i)^{}$'.format(self.repr_def.replace('+', '\+'))

    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.vars_def.keys():
            regex = regex.replace(var, SPEC_REGEX_MATCHING)

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
        'auto-dephosphorylation',
        '$x_ap-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{P} -> $y#$y-{0}'
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
        'deubiquitination',
        '$x_ub-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x# + $y#$y-{ub} -> $x# + $y#$y-{0}'
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
        '$x#$x--0 + $y#$y--0 -> $x,$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'protein-protein-dissociation',
        '$x_ppi-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x,$y#$x--$y -> $x#$x--0 + $y#$y--0'
    ),
    ReactionDef(
        STATE_DEFS,
        'interaction',
        '$x_i+_$y',
        {
            '$x': (Spec, LocusResolution.domain),
            '$y': (Spec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x,$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'dissociation',
        '$x_i-_$y',
        {
            '$x': (Spec, LocusResolution.domain),
            '$y': (Spec, LocusResolution.domain)
        },
        '$x,$y#$x--$y -> $x#$x--0 + $y#$y--0'
    ),
    ReactionDef(
        STATE_DEFS,
        'transcription',
        '$x_trsc_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (GeneSpec, LocusResolution.component)
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
        '$x_ipi+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x--0,$y--0 -> $x#$x--[$y.locus]'
    ),
    ReactionDef(
        STATE_DEFS,
        'intra-protein-dissociation',
        '$x_ipi-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x#$x--[$y.locus] -> $x#$x--0,$y--0'
    ),
    ReactionDef(
        STATE_DEFS,
        'gene-protein-interaction',
        '$x_bind+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (GeneSpec, LocusResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x,$y#$x--$y'
    ),
    ReactionDef(
        STATE_DEFS,
        'gene-protein-dissociation',
        '$x_bind-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (GeneSpec, LocusResolution.domain)
        },
        '$x,$y#$x--$y -> $x#$x--0 + $y#$y--0'
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
    ),
    ReactionDef(
        STATE_DEFS,
        'auto-GuanineNucleotideExchange',
        '$x_agex_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{0} -> $y#$y-{GTP}'
    ),
    ReactionDef(
        STATE_DEFS,
        'auto-GTPHydrolysis',
        '$x_aghy_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{GTP} -> $y#$y-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'GTPase-activation',
        '$x_gap_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x# + $y#$y-{GTP} -> $x# + $y#$y-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'guanine-nucleotide-exchange',
        '$x_gef_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x# + $y#$y-{0} -> $x# + $y#$y-{GTP}'
    ),
    ReactionDef(
        STATE_DEFS,
        'nuclear-export',
        '$x_nexp_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$x# + $y#$y_[(loc)]-{nucleus} -> $x# + $y#$y_[(loc)]-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'nuclear-import',
        '$x_nimp_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$x# + $y#$y_[(loc)]-{0} -> $x# + $y#$y_[(loc)]-{nucleus}'
    ),
    ReactionDef(
        STATE_DEFS,
        'flip-out',
        '$x_FlipOut_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$y#$y_[(direction)]-{0} -> $y#$y_[(direction)]-{out}'
    ),
    ReactionDef(
        STATE_DEFS,
        'flip-in',
        '$x_FlipIn_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$y#$y_[(direction)]-{out} -> $y#$y_[(direction)]-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'facilitated-uptake-extracellular-cytoplasmic',
        '$x_FDExtCyt_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x# + $y#$y_[(loc)]-{0} -> $x# + $y#$y_[(loc)]-{cytosol}'
    ),
    ReactionDef(
        STATE_DEFS,
        'facilitated-uptake-cytoplasmic-extracellular',
        '$x_FDCytExt_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x# + $y#$y_[(loc)]-{cytosol} -> $x# + $y#$y_[(loc)]-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'truncation',
        '$x_CUT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.residue)
        },
        '$x# + $y#$y-{0} -> $x# + $y#$y-{truncated}'
    ),
    ReactionDef(
        STATE_DEFS,
        'satellite',
        '$x_SAT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{0} -> $y#$y-{SAT}'
    ),
    ReactionDef(
        STATE_DEFS,
        'DuplicationPlaque',
        '$x_DUP_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{SAT} -> $y#$y-{DUP}'
    ),
    ReactionDef(
        STATE_DEFS,
        'SPBDuplication',
        '$x_SPB_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{DUP} -> $y#$y-{SPB}'
    ),
    ReactionDef(
        STATE_DEFS,
        'SPBSeparation',
        '$x_SEP_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{SPB} -> $y#$y-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'DNAlicensing',
        '$x_LIC_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{} -> $y#$y-{LIC}'
    ),
    ReactionDef(
        STATE_DEFS,
        'DNAreplicationInitiation',
        '$x_RepInit_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{LIC} -> $y#$y-{RepInit}'
    ),
    ReactionDef(
        STATE_DEFS,
        'DNAreplication',
        '$x_REP_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{RepInit} -> $y#$y-{REP}'
    ),
    ReactionDef(
        STATE_DEFS,
        'DNAsegregation',
        '$x_SEG_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{REP} -> $y#$y-{0}'
    ),
    ReactionDef(
        STATE_DEFS,
        'CYTpolarisation',
        '$x_POL_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{0} -> $y#$y-{POL}'
    ),
    ReactionDef(
        STATE_DEFS,
        'CYTbudding',
        '$x_BUD_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{POL} -> $y#$y-{BUD}'
    ),
    ReactionDef(
        STATE_DEFS,
        'CYTisotropic',
        '$x_ISO_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{BUD} -> $y#$y-{ISO}'
    ),
    ReactionDef(
        STATE_DEFS,
        'CYTokinesis',
        '$x_CYT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y#$y-{ISO} -> $y#$y-{0}'
    )]


class Reaction:
    def __init__(self, definition: ReactionDef, vars: Dict[str, Any]):
        self._definition = definition
        self._vars       = vars

        self.name        = definition.name
        self.terms_lhs   = definition.terms_lhs_from_vars(vars)
        self.terms_rhs   = definition.terms_rhs_from_vars(vars)
        self._repr       = definition.repr_from_vars(vars)

        self._consumed_states    = None
        self._produced_states    = None
        self._synthesised_states = None
        self._degraded_states    = None

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self._repr

    def __eq__(self, other: 'Reaction') -> bool:
        return self._definition == other._definition and self._vars == other._vars

    def invalidate_state_cache(self):
        self._consumed_states    = None
        self._produced_states    = None
        self._synthesised_states = None
        self._degraded_states    = None

    @property
    def components_lhs(self) -> List[Spec]:
        return [spec for term in self.terms_lhs for spec in term.specs]

    @property
    def components_rhs(self) -> List[Spec]:
        return [spec for term in self.terms_rhs for spec in term.specs]

    @property
    def consumed_states(self) -> List[State]:
        if self._consumed_states is None:
            self._consumed_states = []

            for term in self.terms_lhs:
                if all(spec in self.components_rhs for spec in term.specs):
                    self._consumed_states += term.states

        return self._consumed_states

    @property
    def produced_states(self) -> List[State]:
        if self._produced_states  is None:
            self._produced_states = []

            for term in self.terms_rhs:
                if all(spec in self.components_lhs for spec in term.specs):
                    self._produced_states += term.states

        return self._produced_states

    @property
    def degraded_states(self) -> List[State]:
        if self._degraded_states is None:
            self._degraded_states = []

            for term in self.terms_lhs:
                if not any(spec in self.components_rhs for spec in term.specs):
                    self._degraded_states += term.states

        return self._degraded_states

    @property
    def synthesised_states(self) -> List[State]:
        if self._synthesised_states is None:
            self._synthesised_states = []

            for term in self.terms_rhs:
                if not any(spec in self.components_lhs for spec in term.specs):
                    self._synthesised_states += term.states

        return self._synthesised_states

    @property
    def degraded_components(self) -> List[Spec]:
        return [component for component in self.components_lhs if component not in self.components_rhs]

    @property
    def synthesised_components(self) -> List[Spec]:
        return [component for component in self.components_rhs if component not in self.components_lhs]


def matching_reaction_def(repr: str) -> Optional[ReactionDef]:
    return next((reaction_def for reaction_def in REACTION_DEFS if reaction_def.matches_repr(repr)), None)


def reaction_from_str(repr: str, standardize=True) -> Reaction:
    def fixed_spec_types(reaction_def: ReactionDef, vars: Dict[str, Any]) -> Dict[str, Any]:
        keys = vars.keys()
        assert len(list(keys)) == 2

        for key in keys:
            required_type = reaction_def.vars_def[key][0]
            if not isinstance(vars[key], required_type):
                if required_type is GeneSpec:
                    vars[key] = vars[key].to_dna_component_spec()
                elif required_type is MRNASpec:
                    vars[key] = vars[key].to_mrna_component_spec()
                elif required_type is ProteinSpec:
                    vars[key] = vars[key].to_protein_component_spec()
                else:
                    raise NotImplementedError

        return vars

    def fixed_resolutions(reaction_def: ReactionDef, vars: Dict[str, Any]) -> Dict[str, Any]:
        keys = vars.keys()
        assert len(list(keys)) == 2

        for key in keys:
            if reaction_def.vars_def[key][1] < vars[key].resolution:
                raise SyntaxError('In reaction {0}, the specified resolution for variable {1}\n is higher than the required {2}'
                                  .format(repr, str(vars[key]), reaction_def.vars_def[key][1]))

            other = [x for x in keys if x != key][0]
            if not vars[key].has_resolution(reaction_def.vars_def[key][1]):
                if reaction_def.vars_def[key][1] == LocusResolution.domain:
                    vars[key].locus.domain = vars[other].component_name
                elif reaction_def.vars_def[key][1] == LocusResolution.residue:
                    vars[key].locus.residue = vars[other].component_name
                else:
                    raise NotImplementedError

        return vars

    repr = repr.strip()
    reaction_def = matching_reaction_def(repr)

    if not reaction_def:
        raise SyntaxError('Could not match reaction {} with definition'.format(repr))

    vars = reaction_def.vars_from_repr(repr)

    if standardize:
        vars = fixed_spec_types(reaction_def, vars)
        vars = fixed_resolutions(reaction_def, vars)

    reaction_def.validate_vars(vars)
    return Reaction(reaction_def, vars)
