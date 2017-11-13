"""Module containing the classes ReactionTerm, ReactionDef, Reaction. Also contains the default reaction
definitions and the list of bidirectional reactions (that get split into two unidirectional ones, but for
which the contingency is only applied to the forward reaction)."""


from typing import Dict, Any, List, Optional, Callable, TypeVar
import re
from copy import deepcopy

from rxncon.core.spec import Spec, MRNASpec, ProteinSpec, LocusResolution, GeneSpec, spec_from_str
from rxncon.core.state import State, state_from_str, FullyNeutralState


T = TypeVar('T')


SPEC_REGEX_MATCHING = r'([A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
SPEC_REGEX_NON_MATCHING = r'(?:[A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
OUTPUT_REACTION_REGEX = r'^\[[\w-]*\]$'
# These reactions get split into two unidirectional ones, i.e. 'ppi' becomes 'ppi+' and 'ppi-'.
# The contingency then only gets applied on the forward reaction.
BIDIRECTIONAL_REACTIONS = [
    'ppi', 'ipi', 'i', 'bind'
]


class ReactionTerm:
    """ReactionTerm describes a single term in a reaction rule. The `specs` are the components on which the
    `states` live. The list `states` could be empty, e.g. for a kinase. However, the specs should all be
    components, and the states should all be elemental."""
    def __init__(self, specs: List[Spec], states: List[State]) -> None:
        if not all(spec.is_component_spec for spec in specs):
            raise SyntaxError('ReactionTerm contains non-component specs. Specs: {}; States: {}'.format(specs, states))

        if not all(state.is_elemental for state in states if not isinstance(state, FullyNeutralState)):
            raise SyntaxError('ReactionTerm contains non-elemental states. Specs: {}; States: {}'.format(specs, states))

        self.specs, self.states = specs, states

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ReactionTerm):
            return NotImplemented
        return self.specs == other.specs and self.states == other.states

    def __str__(self) -> str:
        return 'ReactionTerm<{0}:{1}>'.format(','.join(str(spec) for spec in self.specs),
                                              ','.join(str(x) for x in self.states))

    def __repr__(self) -> str:
        return str(self)


class ReactionDef:
    """ReactionDef provides meaning for reaction strings such as A_[x]_ppi+_B_[y]. It contains the parsed
    data from the reaction string, and the underlying skeleton rule (`rule_def`)."""
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
        return 'ReactionDef: {0}; name_def: {1}; rule_def: {2} '.format(self.reaction_class, self.name_def,
                                                                        self.rule_def)

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
                'In Reaction {0}: {1} is of type {2}, required to be of type {3}\nVariables: {4}' \
                    .format(self.name_def, var, type(val), self.vars_def[var][0], var_to_val)
            assert val.has_resolution(self.vars_def[var][1]), \
                'In Reaction {0}: {1} is of resolution {2}, required to be of resolution {3}\nVariables: {4}' \
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

        try:
            return ReactionTerm(parse_specs_str(specs_str), parse_states_str(states_str))
        except SyntaxError as e:
            raise SyntaxError('Could not parse reaction {}, {}'.format(term_def, e.msg))

    def _to_base_regex(self) -> str:
        # The (?i) makes the regex case insensitive.
        return r'(?i)^{}$'.format(self.name_def.replace('+', r'\+'))

    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.vars_def.keys():
            regex = regex.replace(var, SPEC_REGEX_MATCHING)

        return regex


# This array will get filled with the ReactionDefs by the function 'initialize_reaction_defs' (below).
# In addition to the DEFAULT_REACTION_DEFS, reactions unique to the system under study can be defined
# in the Excel sheet containing the rxncon system.
REACTION_DEFS = []  # type: List[ReactionDef]


DEFAULT_REACTION_DEFS = [
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
        'protein-gene-interaction',
        '$x_bind+_$y',
        {
            '$x': (ProteinSpec, LocusResolution.domain),
            '$y': (GeneSpec, LocusResolution.domain)
        },
        '$x%#$x%--0 + $y%#$y%--0 -> $x%!$y%#$x%--$y%'
    ),
    ReactionDef(
        'protein-gene-dissociation',
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
        'truncation',
        '$x_CUT_$y',
        {
            '$x': (ProteinSpec, LocusResolution.component),
            '$y': (Spec, LocusResolution.residue)
        },
        '$x%# + $y%#$y%-{0} -> $x%# + $y%#$y%-{truncated}'
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
]


def initialize_reaction_defs(additional_defs: List[Dict[str, str]]=None) -> None:
    """Modifies the module-level variable REACTION_DEFS by adding those reaction definitions inside
    DEFAULT_REACTION_DEFS and additional ones passed into this function."""
    global REACTION_DEFS
    global DEFAULT_REACTION_DEFS
    global BIDIRECTIONAL_REACTIONS

    if not additional_defs:
        additional_defs = []

    type_str_to_spec = {
        'Any': Spec,
        'Protein': ProteinSpec,
        'Gene': GeneSpec,
        'mRNA': MRNASpec,
    }

    res_str_to_spec = {
        'component': LocusResolution.component,
        'domain': LocusResolution.domain,
        'residue': LocusResolution.residue
    }

    parsed_defs = []  # type: List[ReactionDef]

    for additional_def in additional_defs:
        raw_rxn_def = ReactionDef(
            additional_def['!UID:Reaction'],
            '$x_{}_$y'.format(additional_def['!UID:ReactionKey']),
            {
                '$x': (type_str_to_spec[additional_def['!MolTypeX']], res_str_to_spec[additional_def['!ResolutionX']]),
                '$y': (type_str_to_spec[additional_def['!MolTypeY']], res_str_to_spec[additional_def['!ResolutionY']])
            },
            additional_def['!SkeletonRule']
        )

        if additional_def['!BidirectionalVerb'] == 'yes':
            assert additional_def['!UID:ReactionKey'][-1] not in ('+', '-'), \
                'The verb string for a bidirectional reaction cannot end in + or -, error in ReactionDef:' \
                ' {}'.format(additional_def['!UID:Reaction'])

            BIDIRECTIONAL_REACTIONS.append(additional_def['!UID:ReactionKey'])

            forward_def, reverse_def = deepcopy(raw_rxn_def), deepcopy(raw_rxn_def)
            # Relabel the classes.
            forward_def.reaction_class += '-forward'
            reverse_def.reaction_class += '-reverse'

            # Split the verb into forward / reverse.
            forward_def.name_def = '$x_{}+_$y'.format(additional_def['!UID:ReactionKey'])
            reverse_def.name_def = '$x_{}-_$y'.format(additional_def['!UID:ReactionKey'])

            # Reverse the arrow for the reverse rule.
            rev_rule_rhs, rev_rule_lhs = forward_def.rule_def.split('->')
            reverse_def.rule_def = '{} -> {}'.format(rev_rule_lhs, rev_rule_rhs)

            processed_defs = [forward_def, reverse_def]
        else:
            processed_defs = [raw_rxn_def]


        parsed_defs += processed_defs

    REACTION_DEFS = DEFAULT_REACTION_DEFS + parsed_defs


initialize_reaction_defs()


class Reaction:  # pylint: disable=too-many-instance-attributes
    """Generic Reaction class. A Reaction is basically a unidirectional rule with terms on the left hand side
    and terms on the right hand side, which is defined by its ReactionDef. A State is 'consumed' if it is on
    the LHS, but not on the RHS (and its component does not vanish); it is 'produced' if it is on the RHS but
    not on the LHS (and its component does not vanish); it is degraded if its component vanishes going from
    LHS to RHS; it is synthesised if its component appears going from LHS to RHS."""
    def __init__(self, reaction_def: ReactionDef, var_to_val: Dict[str, Any]) -> None:
        self.reaction_class = reaction_def.reaction_class
        self.reaction_def = reaction_def
        self.var_to_val = var_to_val

        self.terms_lhs = reaction_def.terms_lhs_from_vars(var_to_val)
        self.terms_rhs = reaction_def.terms_rhs_from_vars(var_to_val)
        self.name = reaction_def.name_from_vars(var_to_val)

        # This information is cached, and lazily initialised.
        self._consumed_states = None  # type: Optional[List[State]]
        self._produced_states = None  # type: Optional[List[State]]
        self._synthesised_states = None  # type: Optional[List[State]]
        self._degraded_states = None  # type: Optional[List[State]]

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
        self._consumed_states = None
        self._produced_states = None
        self._synthesised_states = None
        self._degraded_states = None

    @property
    def components_lhs(self) -> List[Spec]:
        return [spec for term in self.terms_lhs for spec in term.specs]

    @property
    def components_rhs(self) -> List[Spec]:
        return [spec for term in self.terms_rhs for spec in term.specs]

    @property
    def components(self) -> List[Spec]:
        return list(set(self.components_lhs + self.components_rhs))

    @property
    def components_lhs_structured(self) -> List[Spec]:
        return [spec.with_struct_index(i) for i, spec in enumerate(self.components_lhs)]

    @property
    def components_rhs_structured(self) -> List[Spec]:
        return [spec.with_struct_index(i) for i, spec in enumerate(self.components_rhs)]

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
        if self._produced_states is None:
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
                    # if none of the components from the lhs are in termn_rhs do:
                    self._synthesised_states += term.states


        return self._synthesised_states

    @property
    def degraded_components(self) -> List[Spec]:
        structured_components = [component for component in self.components_lhs_structured if component not in self.components_rhs_structured]
        components = []
        for component in structured_components:
            # destructurize
            component.struct_index = None
            components.append(component)
        return components

    @property
    def synthesised_components(self) -> List[Spec]:
        structured_components= [component for component in self.components_rhs_structured if component not in self.components_lhs_structured]
        components = []
        for component in structured_components:
            # destructurize
            component.struct_index = None
            components.append(component)
        return components

    @property
    def modifier_components(self) -> List[Spec]:
        res = []  # type: List[Spec]
        for term in self.terms_lhs:
            if term in self.terms_rhs:
                res += term.specs

        return res

    @property
    def modifier_states(self) -> List[State]:
        res = []  # type: List[State]
        for term in self.terms_lhs:
            if term in self.terms_rhs:
                res += term.states

        return res


class OutputReaction(Reaction):
    """OutputReactions do not have an underlying rule: they typically describe something like a change
    in a global quantity; therefore they cannot be described by the generic Reaction class."""
    def __init__(self, name: str) -> None:  # pylint: disable=super-init-not-called
        self.name = name
        self.terms_lhs = []
        self.terms_rhs = []

    def __hash__(self) -> int:
        return hash(str(self))

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
    """Given a ReactionDef string, returns the ReactionDef object or None."""
    return next((reaction_def for reaction_def in REACTION_DEFS if reaction_def.matches_name_def(name)),
                None)  # type: ignore


def reaction_from_str(name: str, standardize: bool=True) -> Reaction:
    """Given a Reaction string, returns the Reaction object, or throws a SyntaxError if  not parsable.
    If 'standardize' is True, the parser is more lenient in what it accepts and first calls 'fixed_spec_types'
    and 'fixed_resolutions'."""
    def fixed_spec_types(reaction_def: ReactionDef, var_to_val: Dict[str, Any]) -> Dict[str, Any]:
        """When the Specs appearing in a Reaction are not of the required type, e.g. a ProteinSpec is passed
        instead of a GeneSpec, this this function modifies them into the required type."""
        keys = var_to_val.keys()
        assert len(list(keys)) == 2

        for key in keys:
            required_type = reaction_def.vars_def[key][0]
            if not isinstance(var_to_val[key], required_type) and isinstance(var_to_val[key], ProteinSpec):
                if required_type is GeneSpec:
                    var_to_val[key] = var_to_val[key].to_gene_component_spec()
                elif required_type is MRNASpec:
                    var_to_val[key] = var_to_val[key].to_mrna_component_spec()
                else:
                    raise NotImplementedError

        return var_to_val

    def fixed_resolutions(reaction_def: ReactionDef, var_to_val: Dict[str, Any]) -> Dict[str, Any]:
        """When the Specs appearing in a Reaction are not of the required resolution, e.g. a component ('A') 
        is passed instead of a residue ('A_[(r)]'), this function modifies them into the required resolution."""
        keys = var_to_val.keys()
        assert len(list(keys)) == 2

        for key in keys:
            if reaction_def.vars_def[key][1] < var_to_val[key].resolution:
                raise SyntaxError(
                    'In reaction {0}, the specified resolution for variable {1}\n is higher than the required {2}'
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
