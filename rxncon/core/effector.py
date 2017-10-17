"""Module containing the classes Effector, StateEffector, NaryEffector, AndEffector, OrEffector, NotEffector,
BooleanOperator, BooleanContingencyName, QualSpec, StructEquivalences, TrivialStructEquivalences"""


import re
from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Tuple, Any
from copy import deepcopy
from enum import Enum, unique
import logging

from rxncon.core.spec import spec_from_str, Spec
from rxncon.core.state import State


BOOLEAN_CONTINGENCY_REGEX = '^<.*>$'


LOGGER = logging.getLogger(__name__)


@unique
class BooleanOperator(Enum):
    op_and = 'and'
    op_or = 'or'
    op_not = 'not'
    op_eqv = 'eqv'


class BooleanContingencyName:
    """BooleanContingencyName holds and validates upon construction names such as `<BOOL>`."""
    def __init__(self, name: str) -> None:
        assert re.match(BOOLEAN_CONTINGENCY_REGEX, name), 'BooleanContingencyName {} does not match ' \
                                                          'regex.'.format(name)
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
    """QualSpec holds a Spec that lives in a namespace. The namespace is a list of
    (stringified) boolean contingency names."""
    def __init__(self, namespace: List[str], spec: Spec) -> None:
        self.namespace = namespace
        self.spec = spec
        self._name = '.'.join(namespace + [str(spec)])

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
    spec = spec_from_str(qualified_spec_str.split('.')[-1])

    return QualSpec(namespace, spec)


class StructEquivalences:
    """StructEquivalences holds equivalence classes of QualSpecs. Every element in an equivalence class corresponds
    to the same physical molecule. Elements within the same equivalence class might be labelled with a different
    structure index when living in different namespaces."""
    def __init__(self) -> None:
        self.eq_classes = []  # type: List[List[QualSpec]]

    def __str__(self) -> str:
        s = ''
        for num, eq_class in enumerate(self.eq_classes):
            s = s + 'EqClass {}: {}\n'.format(num, str(eq_class))

        return s

    @property
    def specs(self):
        return list(set(qual_spec.spec.to_component_spec().to_non_struct_spec()
                        for eq_class in self.eq_classes for qual_spec in eq_class))

    def add_equivalence(self, first_qual_spec: QualSpec, second_qual_spec: QualSpec) -> None:
        first_qual_spec, second_qual_spec = first_qual_spec.to_component_qual_spec(), \
                                            second_qual_spec.to_component_qual_spec()

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
            if next((x for x in existing_class if x.to_component_qual_spec() in eq_class), False):
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
                existing_spec = deepcopy(
                    next((x.spec for x in eq_class if x.is_in_root_namespace), None))  # type: ignore
                if existing_spec:
                    existing_spec.locus = deepcopy(qual_spec.spec.locus)
                    return existing_spec

        return None

    def indices_in_root_namespace(self) -> List[int]:
        return [qspec.spec.struct_index for eq_class in self.eq_classes for qspec in eq_class
                if qspec.is_in_root_namespace and qspec.spec.struct_index is not None]


class TrivialStructEquivalences(StructEquivalences):
    """TrivialStructEquivalences describes a situation in which all molecules are inequivalent:
    all namespacing information is ignored."""
    def __init__(self, initial_struct_specs: Dict[Spec, Spec]=None) -> None:
        if not initial_struct_specs:
            self.struct_specs = {}  # type: Dict[Spec, Spec]
        else:
            self.struct_specs = initial_struct_specs

        self.cur_index = 2

    def __str__(self) -> str:
        return 'TrivialStructEquivalences'

    @property
    def specs(self):
        return list(set(spec.to_component_spec().to_non_struct_spec() for spec in self.struct_specs.keys()))

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


class StructCounter:
    """StructCounter is a helper class that generates structure indices. By default it starts counting at 2 since the
    numbers 0 and 1 are reserved for the reactants."""
    def __init__(self, counter_start: Optional[int]=2) -> None:
        if counter_start is None:
            counter_start = 2
        self.value = counter_start

    def increment(self) -> None:
        self.value += 1


class Effector(ABC):
    """Effector is the abstract parent class of the different types of Effector."""
    def __init__(self):
        """All children have the `name` attribute, this code will never be called but exists purely to
        satisy the type checker."""
        self.name = None  # type: Optional[str]

    @property
    @abstractmethod
    def states(self) -> List[State]:
        pass

    @property
    def is_structured(self) -> bool:
        raise NotImplementedError

    @property
    def equivs_specs(self) -> List[Spec]:
        return []

    def collect_global_equivs(self, glob_equivs: StructEquivalences=None,
                              counter: StructCounter=None,
                              cur_namespace: List[str]=None) -> Tuple[StructEquivalences, StructCounter]:
        """Recursively collects all StructEquivalences inside the Effector, without mutating anything."""
        raise NotImplementedError

    def to_global_struct_effector(self, glob_equivs: StructEquivalences=None,
                                  counter: StructCounter=None,
                                  cur_namespace: List[str]=None) -> 'Effector':
        """Returns an effector in which the Specs are recursively re-labelled using a StructEquivalences object."""
        raise NotImplementedError

    @staticmethod
    def _init_to_struct_effector_args(glob_equivs: Optional[StructEquivalences], cur_index: Optional[StructCounter],
                                      cur_namespace: Optional[List[str]]) \
            -> Tuple[StructEquivalences, StructCounter, List[str]]:
        if glob_equivs is None:
            glob_equivs = StructEquivalences()
        if cur_index is None:
            cur_index = StructCounter()
        if cur_namespace is None:
            cur_namespace = []

        return glob_equivs, cur_index, cur_namespace


class StateEffector(Effector):
    """StateEffector holds a State and therefore is a leaf in the Effector."""
    def __init__(self, expr: State, name: Optional[str]=None) -> None:
        self.expr = expr
        self.name = name

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

    def collect_global_equivs(self, glob_equivs: StructEquivalences = None,
                              counter: StructCounter = None,
                              cur_namespace: List[str] = None) -> Tuple[StructEquivalences, StructCounter]:
        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        return glob_equivs, counter

    def to_global_struct_effector(self, glob_equivs: StructEquivalences = None,
                                  counter: StructCounter = None,
                                  cur_namespace: List[str] = None) -> Effector:
        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        state = deepcopy(self.expr)
        updates = {}

        assert not (state.is_homodimer and not state.is_structured), \
            'Please provide structure annotation for homodimer {}'.format(state)

        LOGGER.debug('to_global_struct_effector : Merging {}'.format(str(self)))
        LOGGER.debug('to_global_struct_effector : Equivs {}'.format(glob_equivs))

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
        LOGGER.debug('to_global_struct_effector : Result {}'.format(str(state)))

        return StateEffector(state, name=self.name)

    @staticmethod
    def _generate_index(glob_equivs: StructEquivalences, cur_index: StructCounter) -> int:
        index = cur_index.value
        while index in glob_equivs.indices_in_root_namespace():
            cur_index.increment()
            index = cur_index.value

        return index


class NotEffector(Effector):
    """NotEffector holds another Effector object, and possibly a name, which it can inherit from the
    Boolean Contingency from which it was constructed."""
    def __init__(self, expr: Effector, **kwargs: Optional[str]) -> None:
        try:
            self.name = kwargs['name']
        except KeyError:
            self.name = None
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

    @property
    def equivs_specs(self) -> List[Spec]:
        return self.expr.equivs_specs

    def collect_global_equivs(self, glob_equivs: StructEquivalences = None,
                              counter: StructCounter = None,
                              cur_namespace: List[str] = None) -> Tuple[StructEquivalences, StructCounter]:
        if not self.name:
            raise AssertionError('Cannot collect_global_equivs nameless NotEffectors.')

        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        return self.expr.collect_global_equivs(glob_equivs, counter, cur_namespace + [self.name])

    def to_global_struct_effector(self, glob_equivs: StructEquivalences = None,
                                  counter: StructCounter = None,
                                  cur_namespace: List[str] = None) -> Effector:
        if not self.name:
            raise AssertionError('Cannot to_global_struct_effector nameless NotEffectors.')

        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        return NotEffector(self.expr.to_global_struct_effector(glob_equivs, counter, cur_namespace + [self.name]),
                           name=self.name)


class NaryEffector(Effector, ABC):
    """NaryEffector is an abstract parent class for AndEffector and OrEffector. It can also hold a name which is
    derived from the Boolean Contingency from which it was constructed. It holds StructEquivalences between the
    Specs in its own namespace and in the namespaces of its member Effectors."""
    def __new__(cls, *exprs: Effector, **kwargs):
        assert len(exprs) != 0

        if len(exprs) == 1:
            res = deepcopy(exprs[0])
            try:
                res.name = kwargs['name']
            except KeyError:
                res.name = None

            return res
        else:
            return super().__new__(cls)

    def __init__(self, *exprs: Effector, **kwargs: Any) -> None:
        self.exprs = exprs

        try:
            self.name = kwargs['name']
        except KeyError:
            self.name = None

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

    @property
    def equivs_specs(self) -> List[Spec]:
        return self.equivs.specs + [spec for expr in self.exprs for spec in expr.equivs_specs]

    def collect_global_equivs(self, glob_equivs: StructEquivalences = None,
                              counter: StructCounter = None,
                              cur_namespace: List[str] = None) -> Tuple[StructEquivalences, StructCounter]:
        if not self.name:
            raise AssertionError('Cannot collect_global_equivs nameless NaryEffectors.')

        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        glob_equivs.merge_with(self.equivs, cur_namespace)

        for expr in self.exprs:
            glob_equivs, counter = expr.collect_global_equivs(glob_equivs, counter, cur_namespace + [self.name])

        return glob_equivs, counter

    def to_global_struct_effector(self, glob_equivs: StructEquivalences = None,
                                  counter: StructCounter = None,
                                  cur_namespace: List[str] = None) -> Effector:

        if not self.name:
            raise AssertionError('Cannot to_global_struct_effector nameless NaryEffectors.')

        glob_equivs, counter, cur_namespace = self._init_to_struct_effector_args(glob_equivs, counter, cur_namespace)
        return type(self)(*(x.to_global_struct_effector(
            glob_equivs, counter, cur_namespace + [self.name]) for x in self.exprs), name=self.name)


class AndEffector(NaryEffector):
    """AndEffector describes a logical AND between two or more Effectors. If the AndEffector
    contains only a single effector, we replace the AndEffector's construction with that
    single effector."""
    def __deepcopy__(self, memodict: Dict) -> 'AndEffector':
        """It is required to override __deepcopy__ since we are overriding __new__ with the
        functionality to return an object of a different class."""
        res = AndEffector(*(deepcopy(x) for x in self.exprs))
        res.name = self.name
        res.equivs = deepcopy(self.equivs)

        return res

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
    """OrEffector describes a logical OR between two or more Effectors. If the OrEffector
    contains only a single effector, we replace the OrEffector's construction with that
    single effector."""
    def __deepcopy__(self, memodict: Dict) -> 'OrEffector':
        """It is required to override __deepcopy__ since we are overriding __new__ with the
        functionality to return an object of a different class."""
        res = OrEffector(*(deepcopy(x) for x in self.exprs))
        res.name = self.name
        res.equivs = deepcopy(self.equivs)

        return res

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
