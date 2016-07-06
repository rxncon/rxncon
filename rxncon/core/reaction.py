from typing import Dict, Tuple, List
import re
from typecheck import typecheck

from rxncon.core.specification import Specification, RnaSpecification, ProteinSpecification, SpecificationResolution, \
    DnaSpecification
from rxncon.core.state import State
from rxncon.syntax.rxncon_from_string import state_from_string, specification_from_string


class Reactant:
    def __init__(self, component: Specification, state: State):
        self.component, self.state = component, state

    def __str__(self) -> str:
        return 'Reactant<{0}>:{1}'.format(str(self.component), str(self.state))

    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def __eq__(self, other: 'Reactant'):
        return self.component == other.component and self.state == other.state


class ReactionDefinition:
    ARROW_TWO_HEADS = '<->'
    ARROW_ONE_HEAD = '->'
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'  # substring matched by the group cannot be retrieved after performing a match or referenced later in the pattern.

    def __init__(self, name: str, representation_def: str, variables_def: Dict[str, Tuple], reactants_def: str):
        self.name, self.representation_def, self.variables_def, self.reactants_defs = \
            name, representation_def, variables_def, reactants_def

        self._parse_reactants_def()

    @typecheck
    def __eq__(self, other: 'ReactionDefinition') -> bool:
        return self.name == other.name and self.representation_def == other.representation_def \
            and self.variables_def == other.variables_def and self.reactants_defs == other.reactants_defs

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'ReactionDefinition: {0}; representation_def: {1}; reactants_defs: {2} '\
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
        def parse_state(state_str: str, variables: Dict) -> State:
            for var, val in variables.items():
                var_resolution_str_domain = '{0}.domain'.format(var.replace('$', '\$'))
                if re.search(var_resolution_str_domain, state_str):
                    state_str = re.sub(var_resolution_str_domain, '[{0}]'.format(str(val.domain)), state_str)
                else:
                    state_str = state_str.replace(var, str(val))

            return state_from_string(state_str)

        def parse_component(component_str: str, variables: Dict) -> Specification:
            component_parts = component_str.split('.')

            base_component = variables[component_parts[0]].to_component_specification()

            if len(component_parts) == 1:
                component = base_component
            elif len(component_parts) == 2 and component_parts[1] == 'mRNA':
                component = base_component.to_rna_component_specification()
            elif len(component_parts) == 2 and component_parts[1] == 'gene':
                component = base_component.to_dna_component_specification()
            elif len(component_parts) == 2 and component_parts[1] == 'protein':
                component = base_component.to_protein_component_specification()
            else:
                raise NotImplementedError

            return component

        component_str, state_str = definition.split('#')

        component, state = parse_component(component_str, variables), parse_state(state_str, variables)
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
            val_spec = specification_from_string(val_str)
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
    ReactionDefinition(
        'phosphorylation',
        '$x_p+_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.component),
            '$y': (ProteinSpecification, SpecificationResolution.residue)
        },
        '$x# + $y#$y-{0} -> $x# + $y#$y-{p}'
    ),
    ReactionDefinition(
        'phosphotransfer',
        '$x_pt_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.residue),
            '$y': (ProteinSpecification, SpecificationResolution.residue)
        },
        '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}'
    ),
    ReactionDefinition(
        'protein-protein-interaction',
        '$x_ppi_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.domain),
            '$y': (ProteinSpecification, SpecificationResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y'
    ),

    ReactionDefinition(
        'transcription',
        '$x_trsc_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.component),
            '$y': (DnaSpecification, SpecificationResolution.component)
        },
        '$x# + $y# -> $x# + $y# + $y.mRNA#'
    ),
    ReactionDefinition(
        'translation',
        '$x_trsl_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.component),
            '$y': (RnaSpecification, SpecificationResolution.component)
        },
        '$x# + $y# -> $x# + $y# + $y.protein#'
    ),

    ReactionDefinition(
        'intra-protein-interaction',
        '$x_ipi_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.domain),
            '$y': (ProteinSpecification, SpecificationResolution.domain)
        },
        '$x#$x--0,$y--0 -> $x#$x--$y.domain'
    ),

    ReactionDefinition(
        'gene-protein-interaction',
        '$x_bind_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.domain),
            '$y': (DnaSpecification, SpecificationResolution.domain)
        },
        '$x#$x--0 + $y#$y--0 -> $x#$x--$y + $y#$x--$y'
    )
]


class Reaction:
    def __init__(self, definition: ReactionDefinition, variables):
        self.definition, self.variables = definition, variables

    def __str__(self):
        return self.definition.representation_from_variables(self.variables)

    def __eq__(self, other):
        assert isinstance(other, Reaction)
        return self.definition == other.definition and self.variables == other.variables

    @property
    def reactants_pre(self):
        return self.definition.reactants_pre_from_variables(self.variables)

    @property
    def reactants_post(self):
        return self.definition.reactants_post_from_variables(self.variables)

    @property
    def sources(self):
        sources = []
        for reactant in self.reactants_pre:
            sources += reactant.state

        return list(set(sources))

    @property
    def products(self):
        products = []
        for reactant in self.reactants_post:
            products += reactant.state

        return list(set(products))


def reaction_from_string(reaction_defs: List[ReactionDefinition], representation: str) -> Reaction:
    the_definition = next((reaction_def for reaction_def in reaction_defs
                          if reaction_def.matches_representation(representation)), None)

    assert the_definition, 'Could not match reaction {} with definition'.format(representation)
    variables = the_definition.variables_from_representation(representation)

    return Reaction(the_definition, variables)
