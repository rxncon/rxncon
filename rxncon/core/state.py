from abc import ABCMeta, abstractmethod, abstractproperty
from enum import unique
import typecheck as tc
import re
import typing as tp

from rxncon.util.utils import OrderedEnum
#import rxncon.syntax.string_from_rxncon as sfr
import rxncon.core.specification as spec

from rxncon.syntax.rxncon_from_string import specification_from_string


@unique
class StateModifier(OrderedEnum):
    unmodified = '0'
    phosphor   = 'p'
    ubiquitin  = 'ub'
    guanosintriphosphat = 'gtp'
    truncated  = 'truncated'



class StateDefinition():
    SPEC_REGEX_GROUPED = '([a-zA-Z0-9\/\[\]\(\)_]+)'
    SPEC_REGEX_UNGROUPED = '[a-zA-Z0-9\/\[\]\(\)_]+'

    def __init__(self, name, representation_def, variables_def, superspecification_def):

        self.name, self.representation_def, \
        self.variables_def, self.superspecification_of_def  = name, representation_def, variables_def, superspecification_def

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'state-definition: name={0}, representation_def={1}'.format(self.name, self.representation_def)

    def matches_representation(self, representation):
        return re.match(self._to_matching_regex(), representation)

    def _to_matching_regex(self):
        regex = '{}'.format(self.representation_def)
        for var in self.variables_def.keys():
            regex = regex.replace(var, self.SPEC_REGEX_GROUPED)
        return '^{}$'.format(regex)

    def variables_from_representation(self, representation):
        assert self.matches_representation(representation)
        variables = { }
        for var, var_def in self.variables_def.items():
            var_regex = self.representation_def.replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            if self.variables_def[var] is StateModifier:
                value = state_modifier_from_string(val_str)
            else:
                value = specification_from_string(val_str)

            variables[var] = value

        return variables

    def representation_from_variables(self, variables):
        representation = self.representation_def
        for var, val in variables.items():
            if val is StateModifier:
                representation = representation.replace(var, str(val.value))
            else:
                representation = representation.replace(var, str(val))

        return representation

    def superspec_from_definition(self):
        pass

def state_modifier_from_string(modifier: str):
    return StateModifier(modifier.lower())

# OUTPUT_REGEX = '^\[.+?\]$'
#     INTERACTION_REGEX = '^.+?--.+?$'
#     MODIFIER_REGEX = '.+?-{.+?}'
#     MODIFIER_VALUE_REGEX = '{.+?}'
#     COMPONENT_REGEX = '\w'

STATE_DEFINITION = [
    StateDefinition('interaction-state',
                    '$x--$y',
                    {'$x': spec.Specification,
                     '$y': spec.Specification},
                    ['interaction-state']
                    ),

    StateDefinition('covalent-modification-state',
                    '$x-{$y}',
                    {'$x': spec.Specification,
                     '$y': StateModifier},
                    ['covalent-modification-state']),

    StateDefinition('component-state',
                    '$x',
                    {'$x': spec.Specification},
                    ['interaction-state', 'covalent-modification-state']
                    )

]


class State():
    def __init__(self, definition: StateDefinition, variables):
        self.definition, self.variables = definition, variables

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.definition.representation_from_variables(self.variables)

    def __eq__(self, other):
        return self.definition == other.definition and self.variables == other.variables

    def components(self):
        return [value for value in self.variables.values() if isinstance(value, spec.Specification)]

    def modifier(self):
        return [value for value in self.variables.values() if isinstance(value, StateModifier)]

    def is_superspecification_of(self, other) -> bool:
        superspec = False
        if other.definition.name == self.definition.name or other.definition.name in self.definition.superspecification_of_def:
            for self_component in self.components():
                for other_component in other.components():
                    if self_component.to_component_specification() == other_component.to_component_specification():
                            if self_component.is_superspecification_of(other_component) \
                                    and other_component.is_subspecification_of(self_component):
                                if other.definition.name in self.definition.superspecification_of_def:
                                    superspec = True
                                elif self.modifier() == other.modifier():
                                    superspec = True
                                else:
                                    return False
                            else:
                                return False
        return superspec

    def is_subspecification_of(self, other) -> bool:
        subspec = False
        if other.definition.name == self.definition.name or self.definition.name in other.definition.superspecification_of_def:
            for self_component in self.components():
                for other_component in other.components():
                    if self_component.to_component_specification() == other_component.to_component_specification():
                        if self_component.is_subspecification_of(other_component) \
                                and other_component.is_superspecification_of(self_component):
                            if self.definition.name in other.definition.superspecification_of_def:
                                subspec = True
                            elif self.modifier() == other.modifier():
                                subspec = True
                            else:
                                return False
                        else:
                            return False
        return subspec

        #         if isinstance(other, CovalentModificationState):
        #             return self.modifier == other.modifier and \
        #                     self.substrate.is_subspecification_of(other.substrate)
        #         elif isinstance(other, ComponentState):
        #             return other.is_superspecification_of(self)
        #         return False


def definition_of_state(representation: str):
    the_definition = None
    for definition in STATE_DEFINITION:
        if definition.matches_representation(representation):
            assert not the_definition
            the_definition = definition

    assert the_definition
    return the_definition


def state_from_string(representation: str):
    the_definition = definition_of_state(representation)
    variables = the_definition.variables_from_representation(representation)

    return State(the_definition, variables)



# class State(metaclass=ABCMeta):
#     def __repr__(self) -> str:
#         return str(self)
#
#     @abstractmethod
#     def __hash__(self):
#         pass
#
#     @abstractmethod
#     def __str__(self):
#         pass
#
#     @abstractmethod
#     def is_superspecification_of(self, other) -> bool:
#         pass
#
#     @abstractmethod
#     def is_subspecification_of(self, other) -> bool:
#         pass
#
#     @abstractproperty
#     def components(self) -> List[com.Specification]:
#         pass
#
#
# class ComponentState(State):
#     @tc.typecheck
#     def __init__(self, component: com.Specification):
#         self.component = component
#         self._validate()
#
#     def _validate(self):
#         assert self.component.name is not None
#         assert self.component.domain is None
#         assert self.component.subdomain is None
#         assert self.component.residue is None
#
#     @tc.typecheck
#     def __eq__(self, other: State):
#         return isinstance(other, ComponentState) and self.component == other.component
#
#     def __hash__(self) -> int:
#         return hash("*comp-state-{}".format(self.component))
#
#     def __repr__(self):
#         return str(self)
#
#     def __str__(self):
#         return sfr.string_from_component_state(self)
#
#     @tc.typecheck
#     def is_superspecification_of(self, other: State) -> bool:
#         if isinstance(other, ComponentState):
#             return self.component == other.component
#         elif isinstance(other, InteractionState) or isinstance(other, SelfInteractionState):
#             return self.component.is_superspecification_of(other.first_component) or self.component.is_superspecification_of(other.second_component)
#         elif isinstance(other, CovalentModificationState):
#             return self.component.is_superspecification_of(other.substrate)
#         elif isinstance(other, TranslocationState):
#             return self.component.is_superspecification_of(other.substrate)
#         else:
#             raise NotImplementedError
#
#     @tc.typecheck
#     def is_subspecification_of(self, other: State) -> bool:
#         return isinstance(other, ComponentState) and self.component == other.component
#
#     @property
#     @tc.typecheck
#     def components(self) -> List[com.Specification]:
#         return [self.component]
#
#
# class CovalentModificationState(State):
#     @tc.typecheck
#     def __init__(self, substrate: com.Specification, modifier: StateModifier):
#         self.substrate = substrate
#         self.modifier = modifier
#         self._validate()
#
#     def _validate(self):
#         assert self.substrate.name is not None
# #        assert self.substrate.residue is not None  # the residue has to be defined by the user
#         assert self.modifier.value is not None
#
#     @tc.typecheck
#     def __eq__(self, other: State) -> bool:
#         return isinstance(other, CovalentModificationState) and self.substrate == other.substrate and self.modifier == other.modifier
#
#     def __hash__(self) -> int:
#         return hash('*cov-mod-state-{}-{}*'.format(self.substrate, self.modifier))
#
#     def __str__(self) -> str:
#         return sfr.string_from_covalent_modification_state(self)
#
#     @tc.typecheck
#     def is_superspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, CovalentModificationState) and self.modifier == other.modifier and \
#             self.substrate.is_superspecification_of(other.substrate)
#
#     @tc.typecheck
#     def is_subspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         if isinstance(other, CovalentModificationState):
#             return self.modifier == other.modifier and \
#                     self.substrate.is_subspecification_of(other.substrate)
#         elif isinstance(other, ComponentState):
#             return other.is_superspecification_of(self)
#         return False
#
#     @property
#     @tc.typecheck
#     def components(self) -> List[com.Specification]:
#         return [self.substrate]
#
#
# class InteractionState(State):
#     @tc.typecheck
#     def __init__(self, first_component: com.Specification, second_component: com.Specification):
#         self.first_component = first_component
#         self.second_component = second_component
#         self._validate()
#
#     def _validate(self):
#         assert self.first_component is not None
#         assert self.second_component is not None
#
#     @tc.typecheck
#     def __eq__(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, InteractionState) and self.first_component == other.first_component and \
#             self.second_component == other.second_component
#
#     def __hash__(self) -> int:
#         return hash('*interaction-state-{}-{}*'.format(self.first_component, self.second_component))
#
#     def __str__(self) -> str:
#         return sfr.string_from_inter_protein_interaction_state(self)
#
#     @tc.typecheck
#     def is_superspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, InteractionState) and self.first_component.is_superspecification_of(other.first_component) \
#             and self.second_component.is_superspecification_of(other.second_component)
#
#     @tc.typecheck
#     def is_subspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         if isinstance(other, InteractionState):
#             return self.first_component.is_subspecification_of(other.first_component) \
#             and self.second_component.is_subspecification_of(other.second_component)
#         elif isinstance(other, ComponentState):
#             return other.is_superspecification_of(self) or other.is_superspecification_of(self)
#         return False
#
#     @property
#     @tc.typecheck
#     def components(self) -> List[com.Specification]:
#         return [self.first_component, self.second_component]
#
#
# class SelfInteractionState(State):
#     @tc.typecheck
#     def __init__(self, first_component: com.Specification, second_component: com.Specification):
#         self.first_component = first_component
#         self.second_component = second_component
#         self._validate()
#
#     def _validate(self):
#         assert self.first_component is not None
#         assert self.second_component is not None
#         assert self.first_component.name == self.second_component.name
#
#     @tc.typecheck
#     def __eq__(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, SelfInteractionState) and self.first_component == other.first_component and \
#             self.second_component == other.second_component
#
#     def __hash__(self) -> int:
#         return hash('*self-interaction-state-{}-{}*'.format(self.first_component, self.second_component))
#
#     def __str__(self) -> str:
#         return sfr.string_from_intra_protein_interaction_state(self)
#
#     @tc.typecheck
#     def is_superspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, SelfInteractionState) and self.first_component.is_superspecification_of(other.first_component) \
#             and self.second_component.is_superspecification_of(other.second_component)
#
#     @tc.typecheck
#     def is_subspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         if isinstance(other, SelfInteractionState):
#             return self.first_component.is_subspecification_of(other.first_component) \
#                     and self.second_component.is_subspecification_of(other.second_component)
#         elif isinstance(other, ComponentState):
#             return other.is_superspecification_of(self)
#         return False
#
#     @property
#     @tc.typecheck
#     def components(self) -> List[com.Specification]:
#         return [self.first_component, self.second_component]
#
#
# class TranslocationState(State):
#     @tc.typecheck
#     def __init__(self, substrate: com.Specification, compartment: str):
#         self.substrate = substrate
#         self.compartment = compartment
#         self._validate()
#
#     def _validate(self):
#         assert self.compartment is not None and not ""
#
#     @tc.typecheck
#     def __eq__(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, TranslocationState) and self.substrate == other.substrate and \
#             self.compartment == other.compartment
#
#     def __hash__(self) -> int:
#         return hash('*transloc-state-{}-{}*'.format(self.substrate, self.compartment))
#
#     def __str__(self) -> str:
#         return sfr.string_from_translocation_state(self)
#
#     @tc.typecheck
#     def is_superspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         return isinstance(other, TranslocationState) and self.compartment == other.compartment and \
#             self.substrate.is_superspecification_of(other.substrate)
#
#     @tc.typecheck
#     def is_subspecification_of(self, other: State) -> bool:
#         assert isinstance(other, State)
#         if isinstance(other, TranslocationState):
#             return self.compartment == other.compartment and \
#                     self.substrate.is_subspecification_of(other.substrate)
#         elif isinstance(other, ComponentState):
#             return other.is_superspecification_of(self)
#         return False
#
#     @property
#     @tc.typecheck
#     def components(self) -> List[com.Specification]:
#         return [self.substrate]
#
#
# class GlobalQuantityState(State):
#     @tc.typecheck
#     def __init__(self, name: str):
#         self.name = name
#         self._validate()
#
#     def _validate(self):
#         assert re.match('^\[.+?\]$', self.name)
#
#     def __hash__(self) -> int:
#         return hash('*global-quantity-state-{}*'.format(self.name))
#
#     def __str__(self) -> str:
#         return sfr.string_from_input_state(self)
#
#     def __eq__(self, other):
#         assert isinstance(other, State)
#         return isinstance(other, GlobalQuantityState) and self.components == other.components and self.name == other.name
#
#     def is_produced_by(self, other):
#         """check if a GLobalQuantityState is produced by a certain GlobalQuantityReaction"""
#         pass
#
#     @tc.typecheck
#     def is_superspecification_of(self, other: State) -> bool:
#         return isinstance(other, GlobalQuantityState) and self == other
#
#     @tc.typecheck
#     def is_subspecification_of(self, other: State) -> bool:
#         return isinstance(other, GlobalQuantityState) and self == other
#
#     @property
#     @tc.typecheck
#     def components(self) -> List[com.Specification]:
#         return []
