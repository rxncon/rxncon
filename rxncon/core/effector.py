import re
from abc import ABCMeta, abstractproperty
from typing import List, Optional, Dict, Tuple, Any
from copy import deepcopy
from enum import Enum, unique
import logging

from rxncon.core.spec import spec_from_str, Spec
from rxncon.core.state import State
from rxncon.util.utils import current_function_name

BOOLEAN_CONTINGENCY_REGEX = '^<.*>$'


LOGGER = logging.getLogger(__name__)


@unique
class BooleanOperator(Enum):
    op_and = 'and'
    op_or  = 'or'
    op_not = 'not'
    op_eqv = 'eqv'


class BooleanContingencyName:  # pylint: disable=too-few-public-methods
    def __init__(self, name: str) -> None:
        assert re.match(BOOLEAN_CONTINGENCY_REGEX, name)
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BooleanContingencyName):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.name


class QualSpec:
    def __init__(self, namespace: List[str], spec: Spec) -> None:
        self.namespace = namespace
        self.spec      = spec
        self._name     = '.'.join(namespace + [str(spec)])

    def __str__(self) -> str:
        return self._name

    def __repr__(self) -> str:
        return 'QualSpec<{}>'.format(self._name)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, QualSpec):
            return NotImplemented
        return self.namespace == other.namespace and self.spec == other.spec

    def to_component_qual_spec(self) -> 'QualSpec':
        return QualSpec(self.namespace, self.spec.to_component_spec())

    @property
    def is_in_root_namespace(self) -> bool:
        return not self.namespace

    def with_prepended_namespace(self, extra_namespace: List[str]) -> 'QualSpec':
        new_namespace = deepcopy(extra_namespace) + deepcopy(self.namespace)
        return QualSpec(new_namespace, self.spec)


def qual_spec_from_str(qualified_spec_str: str) -> QualSpec:
    namespace = [x for x in qualified_spec_str.split('.')[:-1]]
    spec      = spec_from_str(qualified_spec_str.split('.')[-1])

    return QualSpec(namespace, spec)


class StructEquivalences:
    def __init__(self) -> None:
        self.eq_classes = []  # type: List[List[QualSpec]]

    def __str__(self) -> str:
        return '\n'.join(str(x) for x in self.eq_classes)

    def add_equivalence(self, first_qual_spec: QualSpec, second_qual_spec: QualSpec) -> None:
        first_qual_spec, second_qual_spec = first_qual_spec.to_component_qual_spec(), second_qual_spec.to_component_qual_spec()

        found_first = None
        found_second = None
        for eq_class in self.eq_classes:
            if first_qual_spec in eq_class:
                if second_qual_spec not in eq_class:
                    eq_class.append(second_qual_spec)
                found_first = eq_class
            elif second_qual_spec in eq_class:
                if first_qual_spec not in eq_class:
                    eq_class.append(first_qual_spec)
                found_second = eq_class

        if found_first and found_second:
            found_first += found_second
            self.eq_classes.remove(found_second)

        if not (found_first or found_second):
            self.eq_classes.append([first_qual_spec, second_qual_spec])

    def add_equivalence_class(self, eq_class: List[QualSpec]) -> None:
        for existing_class in self.eq_classes:
            if next((x for x in existing_class if x.to_component_qual_spec() in eq_class), False):  # type: ignore
                for elem in eq_class:
                    self.add_equivalence(elem.to_component_qual_spec(), existing_class[0])
                return

        self.eq_classes.append([x.to_component_qual_spec() for x in eq_class])

    def merge_with(self, other: 'StructEquivalences', other_base_namespace: List[str]) -> None:
        for other_eq_class in other.eq_classes:
            self.add_equivalence_class([x.with_prepended_namespace(other_base_namespace) for x in other_eq_class])

    def find_unqualified_spec(self, qual_spec: QualSpec) -> Optional[Spec]:
        for eq_class in self.eq_classes:
            if qual_spec.to_component_qual_spec() in eq_class:
                existing_spec = deepcopy(next((x.spec for x in eq_class if x.is_in_root_namespace), None))  # type: ignore
                if existing_spec:
                    existing_spec.locus = deepcopy(qual_spec.spec.locus)
                    return existing_spec

        return None

    def indices_in_root_namespace(self) -> List[int]:
        return [qspec.spec.struct_index for eq_class in self.eq_classes for qspec in eq_class
                if qspec.is_in_root_namespace and qspec.spec.struct_index is not None]


class TrivialStructEquivalences(StructEquivalences):
    def __init__(self, initial_struct_specs: Dict[Spec, Spec]=None) -> None:  # pylint: disable=super-init-not-called
        if not initial_struct_specs:
            self.struct_specs = {}  # type: Dict[Spec, Spec]
        else:
            self.struct_specs = initial_struct_specs

        self.cur_index = 2

    def __str__(self) -> str:
        return 'TrivialStructEquivalences'

    def add_equivalence(self, first_qual_spec: QualSpec, second_qual_spec: QualSpec) -> None:
        pass

    def add_equivalence_class(self, eq_class: List[QualSpec]) -> None:
        pass

    def merge_with(self, other: 'StructEquivalences', other_base_namespace: List[str]) -> None:
        pass

    def find_unqualified_spec(self, qual_spec: QualSpec) -> Optional[Spec]:
        try:
            struct_spec = deepcopy(self.struct_specs[qual_spec.spec.to_component_spec()])
            struct_spec.locus = deepcopy(qual_spec.spec.locus)
            return struct_spec
        except KeyError:
            self.struct_specs[qual_spec.spec.to_component_spec()] = \
                deepcopy(qual_spec.spec.to_component_spec().with_struct_index(self.cur_index))
            self.cur_index += 1
            return self.find_unqualified_spec(qual_spec)

    def indices_in_root_namespace(self) -> List[int]:
        return [x for x in range(self.cur_index)]


class StructCounter:  # pylint: disable=too-few-public-methods
    def __init__(self) -> None:
        self.value = 2

    def increment(self) -> None:
        self.value += 1


class Effector(metaclass=ABCMeta):
    @property
    def name(self) -> Optional[str]:
        try:
            return self._name
        except AttributeError:
            return None

    @name.setter
    def name(self, value: str) -> None:
        self._name = value  # pylint: disable=attribute-defined-outside-init

    @abstractproperty
    def states(self) -> List[State]:
        pass

    @property
    def is_structured(self) -> bool:
        raise NotImplementedError

    def to_merged_struct_effector(self, glob_equivs: StructEquivalences=None,
                                  counter: StructCounter=None,
                                  cur_namespace: List[str]=None) -> 'Effector':
        raise NotImplementedError

    @staticmethod
    def _init_to_struct_effector_args(glob_equivs: Optional[StructEquivalences], cur_index: Optional[StructCounter],
                                      cur_namespace: Optional[List[str]]) -> Tuple[StructEquivalences, StructCounter, List[str]]:
        if glob_equivs is None:
            glob_equivs = StructEquivalences()
        if cur_index is None:
            cur_index = StructCounter()
        if cur_namespace is None:
            cur_namespace = []

        return glob_equivs, cur_index, cur_namespace


class StateEffector(Effector):
    def __init__(self, expr: State) -> None:
        self.expr = expr

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'StateEffector({})'.format(str(self.expr))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Effector):
            return NotImplemented
        return isinstance(other, StateEffector) and self.expr == other.expr and self.name == other.name

    @property
    def states(self) -> List[State]:
        return [deepcopy(self.expr)]

    @property
    def is_structured(self) -> bool:
        return self.expr.is_structured

    def to_merged_struct_effector(self, glob_equivs: StructEquivalences=None,
                                  counter: StructCounter=None,
                                  cur_namespace: List[str]=None) -> Effector:
        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        state = deepcopy(self.expr)
        updates = {}

        assert not (state.is_homodimer and not state.is_structured), 'Please provide structure annotation for homodimer {}'.format(state)

        LOGGER.debug('{} : Merging {}'.format(current_function_name(), str(self)))
        LOGGER.debug('{} : Equivs {}'.format(current_function_name(), glob_equivs))

        for spec in state.specs:
            existing_spec = glob_equivs.find_unqualified_spec(QualSpec(cur_namespace, spec))

            if existing_spec:
                updates[spec] = existing_spec
            else:
                new_spec = deepcopy(spec)
                new_spec.struct_index = self._generate_index(glob_equivs, counter)

                updates[spec] = new_spec
                glob_equivs.add_equivalence(QualSpec([], new_spec), QualSpec(cur_namespace, spec))

        state.update_specs(updates)
        LOGGER.debug('{} : Result {}'.format(current_function_name(), str(state)))

        return StateEffector(state)

    @staticmethod
    def _generate_index(glob_equivs: StructEquivalences, cur_index: StructCounter) -> int:
        index = cur_index.value
        while index in glob_equivs.indices_in_root_namespace():
            cur_index.increment()
            index = cur_index.value

        return index


class NotEffector(Effector):
    def __init__(self, expr: Effector, **kwargs: Optional[str]) -> None:
        try:
            self.name = kwargs['name']
        except KeyError:
            pass
        self.expr = expr

    def __str__(self) -> str:
        return 'NotEffector({})'.format(self.expr)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Effector):
            return NotImplemented
        return isinstance(other, NotEffector) and self.expr == other.expr and self.name == other.name

    @property
    def states(self) -> List[State]:
        return self.expr.states

    @property
    def is_structured(self) -> bool:
        return self.expr.is_structured

    def to_merged_struct_effector(self, glob_equivs: StructEquivalences=None,
                                  counter: StructCounter=None,
                                  cur_namespace: List[str]=None) -> Effector:
        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        return NotEffector(self.expr.to_merged_struct_effector(glob_equivs, counter, cur_namespace), name=self.name)


class NaryEffector(Effector):
    def __init__(self, *exprs: Effector, **kwargs: Any) -> None:
        self.exprs = exprs

        try:
            self.name = kwargs['name']
        except KeyError:
            pass

        try:
            self.equivs = kwargs['equivs']
        except KeyError:
            self.equivs = StructEquivalences()

    @property
    def states(self) -> List[State]:
        return [state for x in self.exprs for state in x.states]

    @property
    def is_structured(self) -> bool:
        return all(x.is_structured for x in self.exprs)

    def to_merged_struct_effector(self, glob_equivs: StructEquivalences=None,
                                  counter: StructCounter=None,
                                  cur_namespace: List[str]=None) -> Effector:
        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        glob_equivs.merge_with(self.equivs, cur_namespace)

        if not self.name:
            raise AssertionError('Cannot merge nameless NaryEffectors.')

        return type(self)(*(x.to_merged_struct_effector(
            glob_equivs, counter, cur_namespace + [self.name]) for x in self.exprs), name=self.name)


class AndEffector(NaryEffector):
    def __str__(self) -> str:
        if self.name:
            return 'AndEffector{0}({1})'.format(self.name, ','.join(str(x) for x in self.exprs))
        else:
            return 'AndEffector({0})'.format(','.join(str(x) for x in self.exprs))

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Effector):
            return NotImplemented
        return isinstance(other, AndEffector) and self.name == other.name and self.exprs == other.exprs


class OrEffector(NaryEffector):
    def __str__(self) -> str:
        if self.name:
            return 'OrEffector{0}({1})'.format(self.name, ','.join(str(x) for x in self.exprs))
        else:
            return 'OrEffector({0})'.format(','.join(str(x) for x in self.exprs))

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Effector):
            return NotImplemented
        return isinstance(other, OrEffector) and self.name == other.name and self.exprs == other.exprs
