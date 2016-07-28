from rxncon.core.specification import Specification, RnaSpecification, ProteinSpecification, SpecificationResolution, DnaSpecification, Domain
from rxncon.core.state import state_from_string, State
from rxncon.syntax.specification_from_string import specification_from_string


from typing import List, Set

import re

class Reactant:
    def __init__(self, component: Specification, states: List[State]):
        self.component, self.states = component, states

    def __str__(self):
        return 'Reactant<{0}>:{1}'.format(str(self.component),
                                          ','.join(str(x) for x in self.states))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.component == other.component and self.states == other.states


def __get_state(state_strs: List[str], variables):
    states = []
    for state_str in state_strs:
        if state_str:
            for var, val in variables.items():
                var_resolution_str_domain = '{0}.domain'.format(var.replace('$', '\$'))
                if re.search(var_resolution_str_domain, state_str):
                    state_str = re.sub(var_resolution_str_domain, '[{0}]'.format(str(val.domain)), state_str)
                else:
                    state_str = state_str.replace(var, str(val))

            states.append(state_from_string(state_str))
    return states

def parse_reactant(definition: str, variables):
    component_str, states_str = definition.split('#')
    component_parts = component_str.split('.')

    if len(component_parts) == 1:
        component = variables[component_parts[0]].to_component_specification()
    elif len(component_parts) == 2 and component_parts[1] == 'mRNA':
        component = variables[component_parts[0]].to_rna_component_specification()
    elif len(component_parts) == 2 and component_parts[1] == 'gene':
        component = variables[component_parts[0]].to_dna_component_specification()
    elif len(component_parts) == 2 and component_parts[1] == 'protein':
        component = variables[component_parts[0]].to_protein_component_specification()
    else:
        raise NotImplementedError

    states = __get_state(states_str.split(','), variables)

    return Reactant(component, states)


class ReactionDefinition:
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'  # substring matched by the group cannot be retrieved after performing a match or referenced later in the pattern.
    def __init__(self, name, representation_def, variables_def, reactants_def):
        self.name, self.representation_def, self.variables_def, self.reactants_defs = \
            name, representation_def, variables_def, reactants_def

        self._parse_reactants_def()

    def __eq__(self, other):
        assert isinstance(other, ReactionDefinition)
        return self.name == other.name and self.representation_def == other.representation_def \
               and self.variables_def == other.variables_def and self.reactants_defs == other.reactants_defs

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'name: {0}; representation_def: {1}; reactants_defs: {2} '.format(self.name, self.representation_def, self.reactants_defs)

    def _parse_reactants_def(self):
        if '<->' in self.reactants_defs:
            arrow = '<->'
        elif '->' in self.reactants_defs:
            arrow = '->'
        else:
            raise AssertionError

        reactants_def_pre_str, reactants_def_post_str = self.reactants_defs.split(arrow)

        self.reactant_defs_pre = [x.strip() for x in reactants_def_pre_str.split('+')]
        self.reactant_defs_post = [x.strip() for x in reactants_def_post_str.split('+')]

    def matches_representation(self, representation):
        return re.match(self._to_matching_regex(), representation)

    def representation_from_variables(self, variables):
        representation = self.representation_def
        for var, val in variables.items():
            representation = representation.replace(var, str(val))

        return representation

    def variables_from_representation(self, representation):
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

    def reactants_pre_from_variables(self, variables):
        return [parse_reactant(x, variables) for x in self.reactant_defs_pre]

    def reactants_post_from_variables(self, variables):
        return [parse_reactant(x, variables) for x in self.reactant_defs_post]

    def _to_base_regex(self):
        return '^{}$'.format(self.representation_def.replace('+', '\+'))

    def _to_matching_regex(self):
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
        'auto-phosphorylation',
        '$x_ap+_$y',
        {
            '$x': (ProteinSpecification, SpecificationResolution.component),
            '$y': (ProteinSpecification, SpecificationResolution.residue)
        },
        '$y#$y-{0} -> $y#$y-{p}'
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
        'protein-protein-interaction',
        '$x_i_$y',
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
    ),

    # ReactionDefinition(
    #     'synthesis',
    #     '$x_syn_$y',
    #     {
    #         '$x': (ProteinSpecification, SpecificationResolution.component),
    #         '$y': (ProteinSpecification, SpecificationResolution.component)
    #     },
    #     '$x# -> $x# + $y#'
    #
    # ),
    # ReactionDefinition(
    #     'degradation',
    #     '$x_deg_$y',
    #     {
    #         '$x': (ProteinSpecification, SpecificationResolution.component),
    #         '$y': (ProteinSpecification, SpecificationResolution.component)
    #     },
    #     '$x# + $y# -> $x#'
    #
    # ),
    #todo: output reaction definition
    # ReactionDefinition(
    #     'output',
    #     '$x',
    #         {
    #             '$x': ()
    #         },
    #     '$x'
    # )
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

def definition_of_state(representation: str):
    the_definition = None

    for definition in REACTION_DEFINITIONS:

        if definition.matches_representation(representation):
            assert not the_definition
            the_definition = definition

    assert the_definition
    return the_definition

def reaction_from_string(representation: str) -> Reaction:
    the_definition = definition_of_state(representation)
    variables = the_definition.variables_from_representation(representation)

    return Reaction(the_definition, variables)

