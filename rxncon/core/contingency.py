"""Module containing the classes Contingency and ContingencyType."""


from enum import unique, Enum
from copy import deepcopy
from collections import defaultdict
import logging
from typing import Any, Callable, Set as TgSet, Dict, List, Optional

from rxncon.core.effector import Effector, StructEquivalences, QualSpec, StateEffector, TrivialStructEquivalences, \
    NotEffector, AndEffector, OrEffector, StructCounter
from rxncon.core.reaction import Reaction, OutputReaction
from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection, Complement, Union, UniversalSet
from rxncon.core.state import State


LOGGER = logging.getLogger(__name__)


@unique
class ContingencyType(Enum):
    """The ContingencyTypes requirement, inhibition are known as `strict` contingencies, whereas the
    positive, negative ContingencyTypes are referred to as quantitative."""
    requirement = '!'
    inhibition = 'x'
    positive = 'k+'
    negative = 'k-'
    no_effect = '0'
    unknown = '?'


class Contingency:
    """Contingency holds the triple `reaction`, `type`, `effector` describing a contingency in a rxncon model.
    Contingency objects are constructed from ContingencyListEntry objects, that live in the module
    rxncon.input.shared.contingency_list."""
    def __init__(self, reaction: Reaction, contingency_type: ContingencyType,
                 effector: Effector, validate_equivs_specs: bool=True) -> None:
        self.reaction, self.contingency_type, self.effector = reaction, contingency_type, effector

        if validate_equivs_specs:
            self.validate_equivs_specs()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Contingency):
            return NotImplemented
        return self.reaction == other.reaction and self.contingency_type == other.contingency_type and \
            self.effector == other.effector

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Contingency({0}, {1}, {2}'.format(str(self.reaction), str(self.contingency_type), str(self.effector))

    def clone(self) -> 'Contingency':
        return deepcopy(self)

    def with_merged_struct_effector(self, equivs: Optional[StructEquivalences]=None,
                                    counter: Optional[StructCounter]=None,
                                    namespace: Optional[List[str]]=None) -> 'Contingency':
        """Returns a Contingency object where the structure information is merged among all Effector objects using
        `equivs`, an object that holds equivalent molecules: different names referring to the same molecule. For
        more details, see the `to_global_struct_effector` and `collect_global_equivs` methods in Effector."""
        structured = self.clone()
        equivs, counter = structured.effector.collect_global_equivs(equivs, counter, namespace)
        structured.effector = structured.effector.to_global_struct_effector(equivs, counter, namespace)
        structured.validate_struct_indices()
        structured.validate_equivs_specs()
        return structured

    def to_structured(self, counter_start: Optional[int]=None) -> 'Contingency':
        """Returns a Contingency object where the structure information is merged among all Effector objects. This
        method first determines the equivalences, after which it calls the `with_merged_struct_effector` method."""
        LOGGER.debug('to_structured: {}'.format(str(self)))

        if isinstance(self.effector, StateEffector) and self.effector.is_structured:
            # A fully structured StateEffector is fine.
            if not self.effector.states[0].is_global and not isinstance(self.reaction, OutputReaction):
                assert any(component in self.reaction.components_lhs_structured for component in
                           self.effector.states[0].components), \
                    "Non-overlapping contingency: {0} does not match structured reaction : {1} (components: {2})" \
                    .format(str(self.effector), str(self.reaction), str(self.reaction.components_lhs_structured))
            return self
        elif isinstance(self.effector, StateEffector) and not self.effector.is_structured:
            # For a non-structured StateEffector, assume the Specs appearing in the Effector
            # match those appearing in the Reaction.
            equivs = StructEquivalences()
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.reaction.components_lhs_structured}
            for spec in self.effector.expr.specs:
                try:
                    equivs.add_equivalence(QualSpec([], struct_components[spec.to_component_spec()]),
                                           QualSpec([str(self.reaction)], spec.to_component_spec()))
                except KeyError:
                    pass

            return self.with_merged_struct_effector(equivs, StructCounter(counter_start), [str(self.reaction)])
        elif self.effector.is_structured:
            # A fully structured Boolean Effector needs to have its structure indices merged.
            return self.with_merged_struct_effector()
        else:
            # For a non-structured Boolean Effector, assume all Specs that could match, actually do match.
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.reaction.components_lhs_structured}
            equivs = TrivialStructEquivalences(struct_components)  # pylint: disable=redefined-variable-type
            return self.with_merged_struct_effector(equivs)

    def to_venn_set(self, k_plus_strict: bool=False, k_minus_strict: bool=False, structured: bool=True,
                    state_wrapper: Callable[[State], Any]=lambda x: x) -> VennSet[Any]:
        """Returns a Venntastic Set object corresponding to the Contingency: requirements are put in a ValueSet,
        inhibitions in a ValueSet within a Complement. If `k_plus_strict` / `k_minus_strict`, then positive and
        negative Contingencies are translated into strict requirements resp. strict inhibitions. If `structured`
        is False, the structure information is discarded. Optionally all States can be wrapped in some other
        class by providing a `state_wrapper`."""
        def parse_effector(eff: Effector) -> VennSet:
            if isinstance(eff, StateEffector):
                if structured:
                    return ValueSet(state_wrapper(eff.expr))
                else:
                    return ValueSet(state_wrapper(eff.expr.to_non_structured()))
            elif isinstance(eff, NotEffector):
                return Complement(parse_effector(eff.expr))
            elif isinstance(eff, OrEffector):
                return Union(*(parse_effector(x) for x in eff.exprs))
            elif isinstance(eff, AndEffector):
                return Intersection(*(parse_effector(x) for x in eff.exprs))
            else:
                raise AssertionError('Unknown Effector {}'.format(str(eff)))

        if k_plus_strict:
            positive = (ContingencyType.requirement, ContingencyType.positive)
        else:
            positive = (ContingencyType.requirement,)  # type: ignore

        if k_minus_strict:
            negative = (ContingencyType.inhibition, ContingencyType.negative)
        else:
            negative = (ContingencyType.inhibition,)  # type: ignore

        if self.contingency_type in positive:
            return parse_effector(self.effector)
        elif self.contingency_type in negative:
            return Complement(parse_effector(self.effector))
        else:
            return UniversalSet()

    def validate_struct_indices(self) -> None:
        """Assert that every index is only used once."""
        specs = [spec for state in self.effector.states for spec in state.specs]
        index_to_specs = defaultdict(set)  # type: Dict[int, TgSet]

        for spec in specs:
            assert spec.struct_index is not None, 'Struct index not assigned in spec {}, ' \
                                                  'contingency {}'.format(spec, self)
            index_to_specs[spec.struct_index].add(spec.to_component_spec())

        assert all(len(x) == 1 for _, x in index_to_specs.items()), 'Structure indices not uniquely assigned in {}'\
            .format(index_to_specs)

    def validate_equivs_specs(self) -> None:
        """Assert that the component Specs appearing in the struct equivalences are appearing either
        in the States of the Effector or in the Reaction."""
        from_equivs = [spec.to_non_struct_spec() for spec in self.effector.equivs_specs]
        from_states = [spec.to_non_struct_spec() for state in self.effector.states for spec in state.components]
        from_reaction = [spec.to_non_struct_spec() for spec in self.reaction.components]

        for spec in from_equivs:
            assert spec in from_states or spec in from_reaction, \
                'Unknown Spec {} appearing in equivalences for reaction {}'.format(spec, self.reaction)






