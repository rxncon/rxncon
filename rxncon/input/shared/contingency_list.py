import re
import logging
from typing import Dict, List, Union, Tuple, Optional
from collections import defaultdict


from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.effector import StateEffector, NotEffector, OrEffector, Effector, AndEffector, \
    BOOLEAN_CONTINGENCY_REGEX, BooleanOperator, BooleanContingencyName, QualSpec, qual_spec_from_str, StructEquivalences
from rxncon.core.reaction import Reaction
from rxncon.core.reaction import reaction_from_str
from rxncon.core.state import state_from_str, State
from rxncon.util.utils import current_function_name

LOGGER = logging.getLogger(__name__)

class ContingencyListEntry:
    def __init__(self, subj: Union[Reaction, BooleanContingencyName],
                 verb: Union[BooleanOperator, ContingencyType],
                 obj: Union[State, BooleanContingencyName, Tuple[QualSpec, QualSpec]]) -> None:
        self.subj = subj
        self.verb = verb
        self.obj  = obj

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ContingencyListEntry):
            return NotImplemented
        return self.subj == other.subj and self.verb == other.verb and self.obj == other.obj

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "ContingencyListEntry<{}, {}, {}>".format(self.subj, self.verb, self.obj)

    @property
    def is_boolean_effector_entry(self) -> bool:
        return self.verb in (BooleanOperator(BooleanOperator.op_and), BooleanOperator(BooleanOperator.op_or),
                             BooleanOperator(BooleanOperator.op_not))

    @property
    def is_boolean_equivalence_entry(self) -> bool:
        return self.verb == BooleanOperator(BooleanOperator.op_eqv)

    @property
    def is_reaction_entry(self) -> bool:
        return isinstance(self.subj, Reaction)


class BooleanContingencyNameWithEquivs(BooleanContingencyName):
    def __init__(self, name: str, equivs: StructEquivalences) -> None:
        super().__init__(name)
        self.equivs = equivs

    def __str__(self) -> str:
        return 'BooleanContingencyNameWithEquivs<{} :: {}>'.format(self.name, self.equivs)


def contingency_list_entry_from_strs(subject_str: str, verb_str: str, object_str: str) -> ContingencyListEntry:
    subject_str, verb_str, object_str = subject_str.strip(), verb_str.lower().strip(), object_str.strip()

    LOGGER.debug('{}: {} / {} / {}'.format(current_function_name(), subject_str, verb_str, object_str))

    subject = None  # type: Optional[Union[Reaction, BooleanContingencyName]]
    verb    = None  # type: Optional[Union[BooleanOperator, ContingencyType]]
    object  = None  # type: Optional[Union[State, BooleanContingencyName, Tuple[QualSpec, QualSpec]]]


    if re.match(BOOLEAN_CONTINGENCY_REGEX, subject_str):
        # subject: Boolean contingency,
        # verb   : Boolean operator,
        # object : State / Boolean contingency / Qual Spec pair.
        subject = BooleanContingencyName(subject_str)
        verb = BooleanOperator(verb_str)
    else:
        # subject: Reaction,
        # verb   : Contingency type,
        # object : State / Boolean contingency.
        subject = reaction_from_str(subject_str)
        verb = ContingencyType(verb_str)

    if re.match(BOOLEAN_CONTINGENCY_REGEX, object_str) and not isinstance(subject, Reaction):
        # subject: Boolean contingency,
        # verb   : Contingency type / Boolean operator,
        # object : Boolean contingency.
        object = BooleanContingencyName(object_str)
    elif re.match(BOOLEAN_CONTINGENCY_REGEX, object_str.split('#')[0]) and isinstance(subject, Reaction):
        # subject: Reaction,
        # verb   : Contingency type,
        # object : Boolean contingency + '#' + reactant equivs.
        name = object_str.split('#')[0]
        equivs_strs = [s.split(',') for s in object_str.split('#')[1:]]
        equivs_dict = {int(i): qual_spec_from_str(qual_spec_str).with_prepended_namespace([name])
                       for i, qual_spec_str in equivs_strs}
        equivs = StructEquivalences()
        for index, spec in enumerate(subject.components_lhs):
            try:
                equivs.add_equivalence(QualSpec([], spec.with_struct_index(index)), equivs_dict[index])
            except KeyError:
                pass

        object = BooleanContingencyNameWithEquivs(name, equivs)
        LOGGER.debug('{} : Created {}'.format(current_function_name(), str(object)))
    elif verb == BooleanOperator.op_eqv:
        strs = [x.strip() for x in object_str.split(',')]
        object = (qual_spec_from_str(strs[0]).with_prepended_namespace([subject.name]),
                  qual_spec_from_str(strs[1]).with_prepended_namespace([subject.name]))
    else:
        object = state_from_str(object_str)

    assert subject is not None
    assert verb is not None
    assert object is not None

    return ContingencyListEntry(subject, verb, object)


def contingencies_from_contingency_list_entries(entries: List[ContingencyListEntry]) -> List[Contingency]:
    contingencies = []

    boolean_entries  = [x for x in entries if x.is_boolean_effector_entry]
    equiv_entries    = [x for x in entries if x.is_boolean_equivalence_entry]
    reaction_entries = [x for x in entries if x.is_reaction_entry]

    effectors        = _create_boolean_contingency_to_effector(boolean_entries)
    equivalences     = _create_boolean_contingency_to_equivalences(equiv_entries)

    while reaction_entries:
        entry = reaction_entries.pop()
        assert isinstance(entry.subj, Reaction)
        contingencies.append(Contingency(entry.subj,
                                         ContingencyType(entry.verb),
                                         _unary_effector_from_boolean_contingency_entry(entry)))

    Effector.dereference = _dereference_boolean_contingency_effectors     # type: ignore
    Effector.contains_booleans = _contains_boolean_contingency_effectors  # type: ignore

    while any(x.effector.contains_booleans() for x in contingencies):     # type: ignore
        for contingency in contingencies:
            contingency.effector.dereference(effectors, equivalences)     # type: ignore

    del Effector.dereference        # type: ignore
    del Effector.contains_booleans  # type: ignore

    return contingencies


class _BooleanContingencyEffector(Effector):
    def __init__(self, expr: BooleanContingencyName, equivs: Optional[StructEquivalences]=None) -> None:
        self.expr = expr
        if not equivs:
            self.equivs = StructEquivalences()
        else:
            self.equivs = equivs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Effector):
            return NotImplemented
        return isinstance(other, _BooleanContingencyEffector) and self.expr == other.expr

    def states(self) -> List[State]:
        return []


def _dereference_boolean_contingency_effectors(self: Effector,
                                               effector_table: Dict[str, Effector],
                                               equivalence_table: Dict[str, List[Tuple[QualSpec, QualSpec]]]) -> None:
    if isinstance(self, _BooleanContingencyEffector):
        LOGGER.debug('{} : {}'.format(current_function_name(), self.expr))
        LOGGER.debug('{} : {}'.format(current_function_name(), effector_table))
        LOGGER.debug('{} : {}'.format(current_function_name(), self.equivs))
        name   = self.expr.name
        equivs = self.equivs
        self.__class__ = effector_table[self.expr.name].__class__
        self.__dict__  = effector_table[self.expr.name].__dict__
        self.name = name
        try:
            self.equivs.merge_with(equivs, [])
            LOGGER.debug('{} : Merged structure information.'.format(current_function_name()))
        except AttributeError:
            self.equivs = equivs
            LOGGER.debug('{} : Initialized structure information.'.format(current_function_name()))
    elif isinstance(self, StateEffector):
        pass
    elif isinstance(self, NotEffector):
        _dereference_boolean_contingency_effectors(self.expr, effector_table, equivalence_table)
    elif isinstance(self, OrEffector) or isinstance(self, AndEffector):
        try:
            assert self.name is not None
            equivs_list = equivalence_table[self.name]

            for equiv in equivs_list:
                self.equivs.add_equivalence(*equiv)

        except KeyError:
            pass

        for expr in self.exprs:
            _dereference_boolean_contingency_effectors(expr, effector_table, equivalence_table)
    else:
        raise AssertionError


def _contains_boolean_contingency_effectors(self: Effector) -> bool:
    if isinstance(self, _BooleanContingencyEffector):
        return True
    elif isinstance(self, StateEffector):
        return False
    elif isinstance(self, NotEffector):
        return _contains_boolean_contingency_effectors(self.expr)
    elif isinstance(self, AndEffector) or isinstance(self, OrEffector):
        return any(_contains_boolean_contingency_effectors(expr) for expr in self.exprs)
    else:
        raise AssertionError


def _create_boolean_contingency_to_effector(boolean_contingencies: List[ContingencyListEntry]) \
        -> Dict[str, Effector]:
    lookup_table = {}  # type: Dict[str, Effector]

    if not boolean_contingencies:
        return lookup_table

    assert all(x.is_boolean_effector_entry for x in boolean_contingencies)

    while boolean_contingencies:
        current_contingency = boolean_contingencies[0]
        assert isinstance(current_contingency.subj, BooleanContingencyName)
        current_contingencies = [x for x in boolean_contingencies if x.subj == current_contingency.subj]
        boolean_contingencies = [x for x in boolean_contingencies if x.subj != current_contingency.subj]

        boolean_operator = BooleanOperator(current_contingency.verb)
        assert all(BooleanOperator(x.verb) == boolean_operator for x in current_contingencies), \
            'Boolean operator inconsistent in contingencies {}'.format(', '.join(str(x) for x in current_contingencies))

        effector_terms = [_unary_effector_from_boolean_contingency_entry(x) for x in current_contingencies]

        if boolean_operator == BooleanOperator.op_and:
            assert len(effector_terms) > 1, 'AND operator {} contains < 2 terms.'.format(' & '.join(str(x) for x in effector_terms))
            effector = AndEffector(*effector_terms)  # type: Effector
        elif boolean_operator == BooleanOperator.op_or:
            assert len(effector_terms) > 1, 'OR operator {} contains < 2 terms.'.format(' & '.join(str(x) for x in effector_terms))
            effector = OrEffector(*effector_terms)
        elif boolean_operator == BooleanOperator.op_not:
            assert len(effector_terms) == 1, 'AND operator {} contains != 1 term.'.format(' & '.join(str(x) for x in effector_terms))
            effector = NotEffector(effector_terms[0])
        else:
            raise AssertionError

        lookup_table[current_contingency.subj.name] = effector

    return lookup_table


def _create_boolean_contingency_to_equivalences(equivalence_contingencies: List[ContingencyListEntry]) \
        -> Dict[str, List[Tuple[QualSpec, QualSpec]]]:
    lookup_table = defaultdict(list)  # type: Dict[str, List[Tuple[QualSpec, QualSpec]]]

    if not equivalence_contingencies:
        return lookup_table

    assert all(x.is_boolean_equivalence_entry for x in equivalence_contingencies)

    for contingency in equivalence_contingencies:
        assert isinstance(contingency.subj, BooleanContingencyName)
        assert isinstance(contingency.obj, tuple)
        lookup_table[contingency.subj.name].append(contingency.obj)

    return lookup_table


def _unary_effector_from_boolean_contingency_entry(entry: ContingencyListEntry) -> Effector:
    if isinstance(entry.obj, State):
        return StateEffector(entry.obj)
    elif isinstance(entry.obj, BooleanContingencyNameWithEquivs):
        return _BooleanContingencyEffector(entry.obj, entry.obj.equivs)
    elif isinstance(entry.obj, BooleanContingencyName):
        return _BooleanContingencyEffector(entry.obj)
    else:
        raise AssertionError
