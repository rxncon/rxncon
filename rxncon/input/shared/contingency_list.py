import re
from functools import reduce
from typing import Dict, List, Union
from typecheck import typecheck

from rxncon.core.reaction import reaction_from_str
from rxncon.util.utils import OrderedEnum
from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.effector import StateEffector, NotEffector, BinaryEffector, OrEffector, Effector, AndEffector
from rxncon.core.reaction import Reaction
from rxncon.core.state import state_from_str, State


BOOLEAN_CONTINGENCY_REGEX = '^<.*>$'


class BooleanOperator(OrderedEnum):
    op_and = 'and'
    op_or  = 'or'
    op_not = 'not'


class BooleanContingencyName:
    @typecheck
    def __init__(self, name: str):
        assert re.match(BOOLEAN_CONTINGENCY_REGEX, name)
        self.name = name

    @typecheck
    def __eq__(self, other: 'BooleanContingencyName') -> bool:
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.name


class ContingencyListEntry:
    @typecheck
    def __init__(self, subject: Union[Reaction, BooleanContingencyName],
                 predicate: Union[BooleanOperator, ContingencyType],
                 agent: Union[State, BooleanContingencyName]):
        self.subject = subject
        self.predicate = predicate
        self.agent = agent

    @typecheck
    def __eq__(self, other: 'ContingencyListEntry') -> bool:
        return self.subject == other.subject and self.predicate == other.predicate and self.agent == other.agent

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "ContingencyListEntry<{}>".format(self.agent)

    @property
    def is_boolean_entry(self) -> bool:
        return isinstance(self.subject, BooleanContingencyName)

    @property
    def is_reaction_entry(self) -> bool:
        return isinstance(self.subject, Reaction)


def contingency_list_entry_from_subject_predicate_agent_strings(subject_str, predicate_str, agent_str) -> ContingencyListEntry:
    predicate_str = predicate_str.lower()

    if re.match(BOOLEAN_CONTINGENCY_REGEX, subject_str):
        subject = BooleanContingencyName(subject_str)
        predicate = BooleanOperator(predicate_str)
    else:
        subject = reaction_from_str(subject_str)
        predicate = ContingencyType(predicate_str)

    if re.match(BOOLEAN_CONTINGENCY_REGEX, agent_str):
        agent = BooleanContingencyName(agent_str)
    else:
        agent = state_from_str(agent_str)

    return ContingencyListEntry(subject, predicate, agent)


@typecheck
def contingencies_from_contingency_list_entries(entries: List[ContingencyListEntry]) -> List[Contingency]:
    contingencies = []

    boolean_entries = [x for x in entries if x.is_boolean_entry]
    reaction_entries = [x for x in entries if x.is_reaction_entry]
    lookup_table = _create_boolean_contingency_lookup_table(boolean_entries)

    while reaction_entries:
        entry = reaction_entries.pop()
        contingencies.append(Contingency(entry.subject,
                                         ContingencyType(entry.predicate),
                                         _unary_effector_from_boolean_contingency_entry(entry)))

    Effector.dereference = _dereference_boolean_contingency_effectors
    Effector.contains_booleans = _contains_boolean_contingency_effectors

    while any([x.effector.contains_booleans() for x in contingencies]):
        for contingency in contingencies:
            contingency.effector.dereference(lookup_table)

    del Effector.dereference
    del Effector.contains_booleans

    return contingencies


class _BooleanContingencyEffector(Effector):
    def __init__(self, expr: BooleanContingencyName):
        self.expr = expr

    @typecheck
    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, _BooleanContingencyEffector) and self.expr == other.expr

    def states(self):
        return [self.expr]


def _dereference_boolean_contingency_effectors(self: Effector, lookup_table: Dict[BooleanContingencyName, Effector]):
    if isinstance(self, _BooleanContingencyEffector):
        name = self.expr.name
        self.__class__ = lookup_table[self.expr].__class__
        self.__dict__ = lookup_table[self.expr].__dict__
        self.name = name
    elif isinstance(self, StateEffector):
        pass
    elif isinstance(self, NotEffector):
        _dereference_boolean_contingency_effectors(self.expr, lookup_table)
    elif isinstance(self, BinaryEffector):
        _dereference_boolean_contingency_effectors(self.left_expr, lookup_table)
        _dereference_boolean_contingency_effectors(self.right_expr, lookup_table)
    else:
        raise AssertionError


def _contains_boolean_contingency_effectors(self: Effector) -> bool:
    if isinstance(self, _BooleanContingencyEffector):
        return True
    elif isinstance(self, StateEffector):
        return False
    elif isinstance(self, NotEffector):
        return _contains_boolean_contingency_effectors(self.expr)
    elif isinstance(self, BinaryEffector):
        return _contains_boolean_contingency_effectors(self.left_expr) or _contains_boolean_contingency_effectors(self.right_expr)
    else:
        raise AssertionError


@typecheck
def _create_boolean_contingency_lookup_table(boolean_contingencies: List[ContingencyListEntry]) -> Dict[BooleanContingencyName, Effector]:
    lookup_table = {}

    if not boolean_contingencies:
        return lookup_table

    assert all([x.is_boolean_entry for x in boolean_contingencies])

    while boolean_contingencies:
        current_contingency = boolean_contingencies[0]
        current_contingencies = [x for x in boolean_contingencies if x.subject == current_contingency.subject]
        boolean_contingencies = [x for x in boolean_contingencies if x.subject != current_contingency.subject]

        boolean_operator = BooleanOperator(current_contingency.predicate)
        assert all([BooleanOperator(x.predicate) == boolean_operator for x in current_contingencies])

        effector_terms = [_unary_effector_from_boolean_contingency_entry(x) for x in current_contingencies]

        if boolean_operator == BooleanOperator.op_and:
            assert len(effector_terms) > 1
            effector = reduce(AndEffector, effector_terms)
        elif boolean_operator == BooleanOperator.op_or:
            assert len(effector_terms) > 1
            effector = reduce(OrEffector, effector_terms)
        elif boolean_operator == BooleanOperator.op_not:
            assert len(effector_terms) == 1
            effector = NotEffector(effector_terms[0])
        else:
            raise AssertionError

        lookup_table[current_contingency.subject] = effector

    return lookup_table


@typecheck
def _unary_effector_from_boolean_contingency_entry(entry: ContingencyListEntry) -> Effector:
    if isinstance(entry.agent, State):
        return StateEffector(entry.agent)
    elif isinstance(entry.agent, BooleanContingencyName):
        return _BooleanContingencyEffector(entry.agent)
    else:
        raise AssertionError
