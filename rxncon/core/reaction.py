from typing import Dict, Any, List, Optional, Callable, TypeVar
import re

from rxncon.core.spec import Spec, MRNASpec, ProteinSpec, LocusResolution, GeneSpec, spec_from_str
from rxncon.core.state import State, state_from_str

T = TypeVar('T')

SPEC_REGEX_MATCHING     = r'([A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
SPEC_REGEX_NON_MATCHING = r'(?:[A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'

OUTPUT_REACTION_REGEX   = r'^\[[\w-]*\]$'

BIDIRECTIONAL_REACTIONS = [
    'ppi', 'ipi', 'i', 'bind'
]


class ReactionTerm:  # pylint:disable=too-few-public-methods
    def __init__(self, specs: List[Spec], states: List[State]) -> None:
        assert all(spec.is_component_spec for spec in specs)
        self.specs, self.states = specs, states

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ReactionTerm):
            return NotImplemented
        return self.states == other.states

    def __str__(self) -> str:
        return 'ReactionTerm<{0}:{1}>'.format(','.join(str(spec) for spec in self.specs), ','.join(str(x) for x in self.states))

    def __repr__(self) -> str:
        return str(self)


class ReactionDef:
    ARROW = '->'

    def __init__(self, reaction_class: str, name_def: str, vars_def: Dict[str, Any], rule_def: str) -> None:
        self.reaction_class, self.name_def, self.vars_def, self.rule_def = reaction_class, name_def, vars_def, rule_def
        self._parse_reactants_def()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ReactionDef):
            return NotImplemented
        return self.reaction_class == other.reaction_class and self.name_def == other.name_def \
            and self.vars_def == other.vars_def and self.rule_def == other.rule_def

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'ReactionDef: {0}; name_def: {1}; rule_def: {2} '.format(self.reaction_class, self.name_def, self.rule_def)

    def matches_name_def(self, name: str) -> bool:
        return True if re.match(self._to_matching_regex(), name) else False

    def name_from_vars(self, var_to_val: Dict[str, Any]) -> str:
        representation = self.name_def
        for var, val in var_to_val.items():
            representation = representation.replace(var, str(val))

        return representation

    def vars_from_name(self, name: str) -> Dict[str, Any]:
        assert self.matches_name_def(name)

        var_to_val = {}
        for var, _ in self.vars_def.items():
            var_regex = self._to_base_regex().replace(var, SPEC_REGEX_MATCHING)
            for other_var in self.vars_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, SPEC_REGEX_NON_MATCHING)

            val_str = re.match(var_regex, name).group(1)
            val = spec_from_str(val_str)

            var_to_val[var] = val

        return var_to_val

    def validate_vars(self, var_to_val: Dict[str, Any]) -> None:
        for var, val in var_to_val.items():
            assert isinstance(val, self.vars_def[var][0]), \
                'In Reaction {0}: {1} is of type {2}, required to be of type {3}\nVariables: {4}'\
                .format(self.name_def, var, type(val), self.vars_def[var][0], var_to_val)
            assert val.has_resolution(self.vars_def[var][1]), \
                'In Reaction {0}: {1} is of resolution {2}, required to be of resolution {3}\nVariables: {4}'\
                .format(self.name_def, var, val.resolution, self.vars_def[var][1], var_to_val)

    def terms_lhs_from_vars(self, var_to_val: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, var_to_val) for x in self.reactant_defs_lhs]

    def terms_rhs_from_vars(self, var_to_val: Dict[str, Any]) -> List[ReactionTerm]:
        return [self._parse_term(x, var_to_val) for x in self.reactant_defs_rhs]

    def _parse_reactants_def(self) -> None:
        assert self.ARROW in self.rule_def

        reactants_def_lhs_str, reactants_def_rhs_str = self.rule_def.split(self.ARROW)

        self.reactant_defs_lhs = [x.strip() for x in reactants_def_lhs_str.split('+')]
        self.reactant_defs_rhs = [x.strip() for x in reactants_def_rhs_str.split('+')]

    def _parse_term(self, term_def: str, var_to_val: Dict[str, Any]) -> ReactionTerm:
        def parse_states_str(states_str: str) -> List[State]:
            if not states_str:
                return []
            try:
                return [parse_vars_str(x, state_from_str) for x in states_str.split('!')]
            except SyntaxError:
                raise SyntaxError('Could not parse states string {}'.format(states_str))

        def parse_specs_str(specs_str: str) -> List[Spec]:
            return [parse_vars_str(x, spec_from_str).to_component_spec() for x in specs_str.split('!')]

        def parse_vars_str(vars_str: str, post_func: Callable[[str], T]) -> T:
            if not vars_str.count('$') == vars_str.count('%'):
                raise SyntaxError('Error in ReactionDef {}'.format(self.reaction_class))

            result_str = vars_str
            for x in re.findall(r'(\$\S*?%)', vars_str):
                var_symbol = x.split('.')[0]

                if var_symbol[-1] == '%':
                    replacement = str(var_to_val[var_symbol[:-1]])
                else:
                    anaphoric_var = var_to_val[var_symbol]  # pylint:disable=unused-variable
                    methods_str = x.split('.', maxsplit=1)[1]
                    assert methods_str[-1] == '%'
                    replacement = str(eval('anaphoric_var.' + methods_str[:-1]))  # pylint:disable=eval-used

                result_str = result_str.replace(x, replacement)

            return post_func(result_str)

        if len(term_def.split('#')) != 2:
            raise SyntaxError('ReactionDef syntax error: each term should be of the form Specs#States,\n'
                              'where Specs is one or more Specs separated by !s and States is zero or\n'
                              'more States separated by !s.\nThe given definition was: {}'.format(term_def))

        specs_str, states_str = term_def.split('#')

        return ReactionTerm(parse_specs_str(specs_str), parse_states_str(states_str))

    def _to_base_regex(self) -> str:
        # The (?i) makes the regex case insensitive.
        return r'(?i)^{}$'.format(self.name_def.replace('+', r'\+'))

    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.vars_def.keys():
            regex = regex.replace(var, SPEC_REGEX_MATCHING)

        return regex


REACTION_DEFS = [
    ReactionDef(
        'phosphorylation',
        '$x_p+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{p}'
    ),
    ReactionDef(
        'dephosphorylation',
        '$x_p-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{p} -> $x%# + $y%#$y%-{0}'
    ),
    ReactionDef(
        'auto-phosphorylation',
        '$x_ap+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{0} -> $y%#$y%-{p}'
    ),
    ReactionDef(
        'auto-dephosphorylation',
        '$x_ap-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{P} -> $y%#$y%-{0}'
    ),
    ReactionDef(
        'ubiquitination',
        '$x_ub+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{ub}'
    ),
    ReactionDef(
        'deubiquitination',
        '$x_ub-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{ub} -> $x%# + $y%#$y%-{0}'
    ),
    ReactionDef(
        'phosphotransfer',
        '$x_pt_$y',
        {
            '$x': (ProteinSpec, LocusResolution.residue),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%#$x%-{p} + $y%#$y%-{0} -> $x%#$x%-{0} + $y%#$y%-{p}'
    ),
    ReactionDef(
        'protein-protein-interaction',
        '$x_ppi+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x%#$x%--0 + $y%#$y%--0 -> $x%!$y%#$x%--$y%'
    ),
    ReactionDef(
        'protein-protein-dissociation',
        '$x_ppi-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x%!$y%#$x%--$y% -> $x%#$x%--0 + $y%#$y%--0'
    ),
    ReactionDef(
        'interaction',
        '$x_i+_$y',
        {
            '$x': (Spec, LocusResolution.domain),
            '$y': (Spec, LocusResolution.domain)
        },
        '$x%#$x%--0 + $y%#$y%--0 -> $x%!$y%#$x%--$y%'
    ),
    ReactionDef(
        'dissociation',
        '$x_i-_$y',
        {
            '$x': (Spec, LocusResolution.domain),
            '$y': (Spec, LocusResolution.domain)
        },
        '$x%!$y%#$x%--$y% -> $x%#$x%--0 + $y%#$y%--0'
    ),
    ReactionDef(
        'transcription',
        '$x_trsc_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (GeneSpec, LocusResolution.component)
        },
        '$x%# + $y%# -> $x%# + $y%# + $y.to_mrna_component_spec()%#0'
    ),
    ReactionDef(
        'translation',
        '$x_trsl_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (MRNASpec, LocusResolution.component)
        },
        '$x%# + $y%# -> $x%# + $y%# + $y.to_protein_component_spec()%#0'
    ),
    ReactionDef(
        'intra-protein-interaction',
        '$x_ipi+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x%#$x%--0!$y%--0 -> $x%#$x%--[$y.locus%]'
    ),
    ReactionDef(
        'intra-protein-dissociation',
        '$x_ipi-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (ProteinSpec, LocusResolution.domain)
        },
        '$x%#$x%--[$y.locus%] -> $x%#$x%--0!$y%--0'
    ),
    ReactionDef(
        'gene-protein-interaction',
        '$x_bind+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (GeneSpec, LocusResolution.domain)
        },
        '$x%#$x%--0 + $y%#$y%--0 -> $x%!$y%#$x%--$y%'
    ),
    ReactionDef(
        'gene-protein-dissociation',
        '$x_bind-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (GeneSpec, LocusResolution.domain)
        },
        '$x%!$y%#$x%--$y% -> $x%#$x%--0 + $y%#$y%--0'
    ),
    ReactionDef(
        'degradation',
        '$x_deg_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x%# + $y%# -> $x%#'
    ),
    ReactionDef(
        'synthesis',
        '$x_syn_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x%# -> $x%# + $y%#0'
    ),
    ReactionDef(
        'auto-GuanineNucleotideExchange',
        '$x_agex_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{0} -> $y%#$y%-{GTP}'
    ),
    ReactionDef(
        'auto-GTPHydrolysis',
        '$x_aghy_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{GTP} -> $y%#$y%-{0}'
    ),
    ReactionDef(
        'GTPase-activation',
        '$x_gap_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{GTP} -> $x%# + $y%#$y%-{0}'
    ),
    ReactionDef(
        'guanine-nucleotide-exchange',
        '$x_gef_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{GTP}'
    ),
    ReactionDef(
        'nuclear-export',
        '$x_nexp_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$x%# + $y%#$y%_[(loc)]-{nucleus} -> $x%# + $y%#$y%_[(loc)]-{0}'
    ),
    ReactionDef(
        'nuclear-import',
        '$x_nimp_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$x%# + $y%#$y%_[(loc)]-{0} -> $x%# + $y%#$y%_[(loc)]-{nucleus}'
    ),
    ReactionDef(
        'flip-out',
        '$x_FlipOut_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$y%#$y%_[(direction)]-{0} -> $y%#$y%_[(direction)]-{out}'
    ),
    ReactionDef(
        'flip-in',
        '$x_FlipIn_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.component)
        },
        '$y%#$y%_[(direction)]-{out} -> $y%#$y%_[(direction)]-{0}'
    ),
    ReactionDef(
        'facilitated-uptake-extracellular-cytoplasmic',
        '$x_FDExtCyt_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x%# + $y%#$y%_[(loc)]-{0} -> $x%# + $y%#$y%_[(loc)]-{cytosol}'
    ),
    ReactionDef(
        'facilitated-uptake-cytoplasmic-extracellular',
        '$x_FDCytExt_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.component)
        },
        '$x%# + $y%#$y%_[(loc)]-{cytosol} -> $x%# + $y%#$y%_[(loc)]-{0}'
    ),
    ReactionDef(
        'truncation',
        '$x_CUT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{truncated}'
    ),
    ReactionDef(
        'satellite',
        '$x_SAT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{0} -> $y%#$y%-{SAT}'
    ),
    ReactionDef(
        'DuplicationPlaque',
        '$x_DUP_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{SAT} -> $y%#$y%-{DUP}'
    ),
    ReactionDef(
        'SPBDuplication',
        '$x_SPB_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{DUP} -> $y%#$y%-{SPB}'
    ),
    ReactionDef(
        'SPBSeparation',
        '$x_SEP_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{SPB} -> $y%#$y%-{0}'
    ),
    ReactionDef(
        'DNAlicensing',
        '$x_LIC_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{0} -> $y%#$y%-{LIC}'
    ),
    ReactionDef(
        'DNAreplicationInitiation',
        '$x_RepInit_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{LIC} -> $y%#$y%-{RepInit}'
    ),
    ReactionDef(
        'DNAreplication',
        '$x_REP_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{RepInit} -> $y%#$y%-{REP}'
    ),
    ReactionDef(
        'DNAsegregation',
        '$x_SEG_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{REP} -> $y%#$y%-{0}'
    ),
    ReactionDef(
        'CYTpolarisation',
        '$x_POL_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{0} -> $y%#$y%-{POL}'
    ),
    ReactionDef(
        'CYTbudding',
        '$x_BUD_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{POL} -> $y%#$y%-{BUD}'
    ),
    ReactionDef(
        'CYTisotropic',
        '$x_ISO_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{BUD} -> $y%#$y%-{ISO}'
    ),
    ReactionDef(
        'CYTokinesis',
        '$x_CYT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$y%#$y%-{ISO} -> $y%#$y%-{0}'
    ),
    ReactionDef(
        'acetylation',
        '$x_Ac+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{Ac}'
    ),
    ReactionDef(
        'Myristoylation',
        '$x_Myr+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{Myr}'
    ),
    ReactionDef(
        'deacetylation',
        '$x_Ac-_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{Ac} -> $x%# + $y%#$y%-{0}'
    ),
    ReactionDef(
        'demyrstoylation',
        '$x%_Myr-_$y%',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (ProteinSpec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{Myr} -> $x%# + $y%#$y%-{0}'
    ),
    ReactionDef(
        'pro-cat-translation',
        '$x_trslprocat_$y',
        {
            '$x': (Spec, LocusResolution.component),
            '$y': (MRNASpec, LocusResolution.component)
        },
        '$x%# + $y%# -> $x%# + $y%# + $y.to_protein_component_spec().with_name_suffix(\'PRO\')%!$y.to_protein_component_spec().with_name_suffix(\'CAT\')%#'                       # pylint: disable=line-too-long
        '$y.to_protein_component_spec().with_name_suffix(\'PRO\').with_domain(\'PROCAT\')%--$y.to_protein_component_spec().with_name_suffix(\'CAT\').with_domain(\'CATPRO\')%!0'  # pylint: disable=line-too-long
    )
]


class Reaction:  # pylint: disable=too-many-instance-attributes
    def __init__(self, reaction_def: ReactionDef, var_to_val: Dict[str, Any]) -> None:
        self.reaction_class = reaction_def.reaction_class
        self.reaction_def   = reaction_def
        self.var_to_val     = var_to_val

        self.terms_lhs = reaction_def.terms_lhs_from_vars(var_to_val)
        self.terms_rhs = reaction_def.terms_rhs_from_vars(var_to_val)
        self.name      = reaction_def.name_from_vars(var_to_val)

        self._consumed_states    = None  # type: Optional[List[State]]
        self._produced_states    = None  # type: Optional[List[State]]
        self._synthesised_states = None  # type: Optional[List[State]]
        self._degraded_states    = None  # type: Optional[List[State]]

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Reaction):
            return NotImplemented
        return self.reaction_def == other.reaction_def and self.var_to_val == other.var_to_val

    def __getitem__(self, item: str) -> Any:
        return self.var_to_val[item]

    def invalidate_state_cache(self) -> None:
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
    def components_lhs_structured(self) -> List[Spec]:
        return [spec.with_struct_index(i) for i, spec in enumerate(self.components_lhs)]

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


class OutputReaction(Reaction):
    def __init__(self, name: str) -> None:  # pylint: disable=super-init-not-called
        self.name      = name
        self.terms_lhs = []
        self.terms_rhs = []

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Reaction):
            return NotImplemented
        return isinstance(other, OutputReaction) and self.name == other.name

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return 'OutputReaction<{}>'.format(self.name)

    @property
    def components_rhs(self) -> List[Spec]:
        return []

    @property
    def consumed_states(self) -> List[State]:
        return []

    @property
    def degraded_states(self) -> List[State]:
        return []

    @property
    def components_lhs(self) -> List[Spec]:
        return []

    @property
    def synthesised_states(self) -> List[State]:
        return []

    @property
    def degraded_components(self) -> List[Spec]:
        return []

    @property
    def synthesised_components(self) -> List[Spec]:
        return []

    @property
    def produced_states(self) -> List[State]:
        return []

    def invalidate_state_cache(self) -> None:
        pass


def matching_reaction_def(name: str) -> Optional[ReactionDef]:
    return next((reaction_def for reaction_def in REACTION_DEFS if reaction_def.matches_name_def(name)), None)  # type: ignore


def reaction_from_str(name: str, standardize: bool=True) -> Reaction:
    def fixed_spec_types(reaction_def: ReactionDef, var_to_val: Dict[str, Any]) -> Dict[str, Any]:
        keys = var_to_val.keys()
        assert len(list(keys)) == 2

        for key in keys:
            required_type = reaction_def.vars_def[key][0]
            if not isinstance(var_to_val[key], required_type) and isinstance(var_to_val[key], ProteinSpec):
                if required_type is GeneSpec:
                    var_to_val[key] = var_to_val[key].to_dna_component_spec()
                elif required_type is MRNASpec:
                    var_to_val[key] = var_to_val[key].to_mrna_component_spec()
                else:
                    raise NotImplementedError

        return var_to_val

    def fixed_resolutions(reaction_def: ReactionDef, var_to_val: Dict[str, Any]) -> Dict[str, Any]:
        keys = var_to_val.keys()
        assert len(list(keys)) == 2

        for key in keys:
            if reaction_def.vars_def[key][1] < var_to_val[key].resolution:
                raise SyntaxError('In reaction {0}, the specified resolution for variable {1}\n is higher than the required {2}'
                                  .format(name, str(var_to_val[key]), reaction_def.vars_def[key][1]))

            other = [x for x in keys if x != key][0]
            if not var_to_val[key].has_resolution(reaction_def.vars_def[key][1]):
                if reaction_def.vars_def[key][1] == LocusResolution.domain:
                    var_to_val[key].locus.domain = var_to_val[other].name
                elif reaction_def.vars_def[key][1] == LocusResolution.residue:
                    var_to_val[key].locus.residue = var_to_val[other].name
                else:
                    raise NotImplementedError

        return var_to_val

    name = name.strip()

    if re.match(OUTPUT_REACTION_REGEX, name):
        return OutputReaction(name)

    reaction_def = matching_reaction_def(name)

    if not reaction_def:
        raise SyntaxError('Could not match reaction {} with definition'.format(name))

    var_to_val = reaction_def.vars_from_name(name)

    if standardize:
        var_to_val = fixed_spec_types(reaction_def, var_to_val)
        var_to_val = fixed_resolutions(reaction_def, var_to_val)

    reaction_def.validate_vars(var_to_val)
    return Reaction(reaction_def, var_to_val)
