"""Module containing ContingencyListEntry, BooleanContingencyNameWithEquivs and the constructor functions
contingency_list_entry_from_strs and contingencies_from_contingency_list_entries. ContingencyListEntries
are the triples (target, type, effector) that appear in the tabular data format used to input rxncon.
In this code, nested boolean expressions are `dereferenced`, such that at the end of the day we always have
a Reaction as a target, doing away with the boolean contingencies that appear as targets.
So-called `equivalences` are structured specs (e.g. A@1) that might have different numerical indices since
they arised from different boolean contingencies, but actually refer to the same molecule."""

import re
import logging
from typing import Dict, List, Union, Tuple, Optional

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.effector import StateEffector, NotEffector, OrEffector, Effector, AndEffector, \
    BOOLEAN_CONTINGENCY_REGEX, BooleanOperator, BooleanContingencyName, QualSpec, qual_spec_from_str, StructEquivalences
from rxncon.core.reaction import Reaction, reaction_from_str
from rxncon.core.spec import Spec
from rxncon.core.state import state_from_str, State


LOGGER = logging.getLogger(__name__)


class ContingencyListEntry:
    """ContingencyListEntry holds the triple (subject, verb, object) or (target, type, effector) appearing
    in the tabular representation of a rxncon system."""
    def __init__(self, subj: Union[Reaction, BooleanContingencyName],
                 verb: Union[BooleanOperator, ContingencyType],
                 obj: Union[State, BooleanContingencyName]) -> None:
        self.subj = subj
        self.verb = verb
        self.obj = obj

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ContingencyListEntry):
            return NotImplemented
        return self.subj == other.subj and self.verb == other.verb and self.obj == other.obj

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "ContingencyListEntry<{}, {}, {}>".format(self.subj, self.verb, self.obj)

    @property
    def is_boolean_entry(self) -> bool:
        return isinstance(self.subj, BooleanContingencyName)

    @property
    def is_reaction_entry(self) -> bool:
        return isinstance(self.subj, Reaction)


class BooleanContingencyNameWithEquivs(BooleanContingencyName):
    """A BooleanContingencyName, carrying with it a set of structure Equivalences."""
    def __init__(self, name: str, equivs: StructEquivalences) -> None:
        super().__init__(name)
        self.equivs = equivs

    def __str__(self) -> str:
        return 'BooleanContingencyNameWithEquivs<{} :: {}>'.format(self.name, self.equivs)


def contingency_list_entry_from_strs(subject_str: str, verb_str: Union[str, float],
                                     object_str: str) -> ContingencyListEntry:
    """Parses a ContingencyListEntry from a triple obtained from a tabular representation of the rxncon system."""
    # The excel parser returns a value of 0, which is used to denote a neutral contingency as a float object.
    if isinstance(verb_str, float):
        verb_str = str(int(verb_str))
        assert verb_str == '0', 'Unrecognized contingency {}'.format(verb_str)

    def _add_equivs(equivs: StructEquivalences, equivs_strs: List[List[str]], name: str) -> StructEquivalences:
        for target_qual_spec_str, source_qual_spec_str in equivs_strs:
            lhs_qual_spec = qual_spec_from_str(target_qual_spec_str)
            rhs_qual_spec = qual_spec_from_str(source_qual_spec_str).with_prepended_namespace([name])
            equivs.add_equivalence(lhs_qual_spec, rhs_qual_spec)
        return equivs

    subject_str, verb_str, object_str = subject_str.strip(), verb_str.lower().strip(), object_str.strip()

    LOGGER.debug('contingency_list_entry_from_strs: {} / {} / {}'.format(subject_str, verb_str, object_str))

    subject = None  # type: Optional[Union[Reaction, BooleanContingencyName]]
    verb = None  # type: Optional[Union[BooleanOperator, ContingencyType]]
    object = None  # type: Optional[Union[State, BooleanContingencyName, Tuple[QualSpec, QualSpec]]]

    if re.match(BOOLEAN_CONTINGENCY_REGEX, subject_str):
        # subject: Boolean contingency,
        # verb   : Boolean operator,
        # object : State / Boolean contingency
        subject = BooleanContingencyName(subject_str)
        verb = BooleanOperator(verb_str)
    else:
        # subject: Reaction,
        # verb   : Contingency type,
        # object : State / Boolean contingency.
        subject = reaction_from_str(subject_str)
        verb = ContingencyType(verb_str)

    if re.match(BOOLEAN_CONTINGENCY_REGEX, object_str) and '#' not in object_str:
        # subject: Boolean contingency, Reaction
        # verb   : Contingency type / Boolean operator,
        # object : Boolean contingency.
        object = BooleanContingencyName(object_str)
    elif re.match(BOOLEAN_CONTINGENCY_REGEX, object_str.split('#')[0]):
        # subject: Reaction / Boolean contingency
        # verb   : Contingency type / Boolean operator
        # object : Boolean contingency + '#' + reactant equivs / Boolean equivs.
        name = object_str.split('#')[0]
        equivs_strs = [s.split('=') for s in object_str.split('#')[1:]]
        equivs = StructEquivalences()

        _add_equivs(equivs, equivs_strs, name)

        object = BooleanContingencyNameWithEquivs(name, equivs)
        LOGGER.debug('contingency_list_entry_from_strs : Created {}'.format(str(object)))
    else:
        object = state_from_str(object_str)

    assert subject is not None, 'Could not parse subject in {} {} {}'.format(subject_str, verb_str, object_str)
    assert verb is not None, 'Could not parse verb in {} {} {}'.format(subject_str, verb_str, object_str)
    assert object is not None, 'Could not parse object in {} {} {}'.format(subject_str, verb_str, object_str)

    return ContingencyListEntry(subject, verb, object)


def contingencies_from_contingency_list_entries(entries: List[ContingencyListEntry]) -> List[Contingency]:
    """Constructs contingencies from the tabular entries. Will dereference pointers to boolean contingencies
    into nested expressions."""
    contingencies = []

    boolean_entries = [x for x in entries if x.is_boolean_entry]
    reaction_entries = [x for x in entries if x.is_reaction_entry]

    while reaction_entries:
        entry = reaction_entries.pop()
        assert isinstance(entry.subj, Reaction)
        contingencies.append(Contingency(entry.subj,
                                         ContingencyType(entry.verb),
                                         _unary_effector_from_boolean_contingency_entry(entry),
                                         validate_equivs_specs=False))

    Effector.dereference = _dereference_boolean_contingency_effectors  # type: ignore
    Effector.contains_booleans = _contains_boolean_contingency_effectors  # type: ignore
    effectors = _create_boolean_contingency_to_effector(boolean_entries)

    while any(x.effector.contains_booleans() for x in contingencies):  # type: ignore
        for contingency in contingencies:
            contingency.effector.dereference(effectors)  # type: ignore

    del Effector.dereference  # type: ignore
    del Effector.contains_booleans  # type: ignore

    for contingency in contingencies:
        contingency.validate_equivs_specs()

    return contingencies


class _BooleanContingencyEffector(Effector):
    def __init__(self, expr: BooleanContingencyName, equivs: Optional[StructEquivalences] = None) -> None:
        self.expr = expr
        if not equivs:
            self.equivs = StructEquivalences()
        else:
            self.equivs = equivs

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Effector):
            return NotImplemented
        return isinstance(other, _BooleanContingencyEffector) and self.expr == other.expr

    @property
    def equivs_specs(self) -> List[Spec]:
        return self.equivs.specs

    @property
    def states(self) -> List[State]:
        return []


def _dereference_boolean_contingency_effectors(self: Effector, effector_table: Dict[str, Effector]) -> None:
    if isinstance(self, _BooleanContingencyEffector):
        LOGGER.debug('_dereference_boolean_contingency_effectors : Expr: {}'.format(self.expr))
        LOGGER.debug('_dereference_boolean_contingency_effectors : EffTable: {}'.format(effector_table))
        LOGGER.debug('_dereference_boolean_contingency_effectors : Equivs: {}'.format(self.equivs))
        name = self.expr.name
        equivs = self.equivs
        self.__class__ = effector_table[self.expr.name].__class__
        self.__dict__ = effector_table[self.expr.name].__dict__
        self.name = name
        try:
            self.equivs.merge_with(equivs, [])
            LOGGER.debug('_dereference_boolean_contingency_effectors : Merged structure information.')
        except AttributeError:
            self.equivs = equivs
            LOGGER.debug('_dereference_boolean_contingency_effectors : Initialized structure information.')
    elif isinstance(self, StateEffector):
        pass
    elif isinstance(self, NotEffector):
        _dereference_boolean_contingency_effectors(self.expr, effector_table)
    elif isinstance(self, OrEffector) or isinstance(self, AndEffector):
        for expr in self.exprs:
            _dereference_boolean_contingency_effectors(expr, effector_table)
    else:
        raise AssertionError('Unknown effector type {} in _dereference_boolean_contingency_effectors'.format(self))


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

    assert all(x.is_boolean_entry for x in boolean_contingencies)

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
            # assert len(effector_terms) > 1, 'AND operator {} contains < 2 terms.'.format(
            #     ' & '.join(str(x) for x in effector_terms))
            effector = AndEffector(*effector_terms)  # type: Effector
        elif boolean_operator == BooleanOperator.op_or:
            # assert len(effector_terms) > 1, 'OR operator {} contains < 2 terms.'.format(
            #     ' & '.join(str(x) for x in effector_terms))
            effector = OrEffector(*effector_terms)
        elif boolean_operator == BooleanOperator.op_not:
            assert len(effector_terms) == 1, 'NOT operator {} contains != 1 term.'.format(
                ' & '.join(str(x) for x in effector_terms))
            effector = NotEffector(effector_terms[0])
        else:
            raise AssertionError

        lookup_table[current_contingency.subj.name] = effector

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
