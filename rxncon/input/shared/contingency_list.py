import re
from enum import Enum
from functools import reduce
from typing import Dict, List

import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.syntax.rxncon_from_string as fst

BOOLEAN_CONTINGENCY_REGEX = '^<.*>$'


class BooleanOperator(Enum):
    op_and = 'and'
    op_or  = 'or'
    op_not = 'not'


class ContingencyListEntry:
    # @todo Names subject, predicate, agent need to be standardized. Maybe validate.
    def __init__(self, subject, predicate, agent):
        self.subject = subject
        self.predicate = predicate
        self.agent = agent

    def __eq__(self, other: 'ContingencyListEntry') -> bool:
        assert isinstance(other, ContingencyListEntry)
        return self.subject == other.subject and self.predicate == other.predicate and self.agent == other.agent

    @property
    def is_boolean_entry(self) -> bool:
        return isinstance(self.subject, BooleanContingencyName)

    @property
    def is_reaction_entry(self) -> bool:
        return isinstance(self.subject, rxn.Reaction)


class BooleanContingencyName:
    def __init__(self, name: str):
        assert re.match(BOOLEAN_CONTINGENCY_REGEX, name)
        self.name = name

    def __eq__(self, other: 'BooleanContingencyName') -> bool:
        assert isinstance(other, BooleanContingencyName)
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)


def contingency_list_entry_from_subject_predicate_agent_strings(subject_string, predicate_string, agent_string) -> ContingencyListEntry:
    # @todo Input/output states, raise Exception when sub-pred-ag is inconsistent triple.
    predicate_string = predicate_string.lower()

    if re.match(BOOLEAN_CONTINGENCY_REGEX, subject_string):
        subject = BooleanContingencyName(subject_string)
        predicate = BooleanOperator(predicate_string)

    else:
        subject = fst.reaction_from_string(subject_string)
        predicate = con.ContingencyType(predicate_string)

    if re.match(BOOLEAN_CONTINGENCY_REGEX, agent_string):
        agent = BooleanContingencyName(agent_string)

    else:
        agent = fst.state_from_string(agent_string)

    return ContingencyListEntry(subject, predicate, agent)


def contingencies_from_contingency_list_entries(entries: List[ContingencyListEntry]) -> List[con.Contingency]:
    # @todo Explain this algorithm (the monkey patching etc.).
    contingencies = []

    boolean_entries = [x for x in entries if x.is_boolean_entry]
    reaction_entries = [x for x in entries if x.is_reaction_entry]
    lookup_table = _create_boolean_contingency_lookup_table(boolean_entries)

    while reaction_entries:
        entry = reaction_entries.pop()
        contingencies.append(con.Contingency(entry.subject,
                                             con.ContingencyType(entry.predicate),
                                             _unary_effector_from_boolean_contingency_entry(entry)))

    eff.Effector.dereference = _dereference_boolean_contingency_effectors
    eff.Effector.contains_booleans = _contains_boolean_contingency_effectors

    while any([x.effector.contains_booleans() for x in contingencies]):
        for contingency in contingencies:
            contingency.effector.dereference(lookup_table)

    del eff.Effector.dereference
    del eff.Effector.contains_booleans

    return contingencies


class _BooleanContingencyEffector(eff.Effector):
    def __init__(self, expr: BooleanContingencyName):
        self.expr = expr

    def __eq__(self, other: eff.Effector) -> bool:
        assert isinstance(other, eff.Effector)

        if isinstance(other, _BooleanContingencyEffector):
            return self.expr == other.expr

        else:
            return False


def _dereference_boolean_contingency_effectors(self: eff.Effector, lookup_table: Dict[BooleanContingencyName, eff.Effector]):
    if isinstance(self, _BooleanContingencyEffector):
        name = self.expr.name
        self.__class__ = lookup_table[self.expr].__class__
        self.__dict__ = lookup_table[self.expr].__dict__
        self.name = name

    elif isinstance(self, eff.StateEffector):
        pass

    elif isinstance(self, eff.NotEffector):
        _dereference_boolean_contingency_effectors(self.expr, lookup_table)

    elif isinstance(self, eff.BinaryEffector):
        _dereference_boolean_contingency_effectors(self.left_expr, lookup_table)
        _dereference_boolean_contingency_effectors(self.right_expr, lookup_table)

    else:
        raise AssertionError


def _contains_boolean_contingency_effectors(self: eff.Effector) -> bool:
    if isinstance(self, _BooleanContingencyEffector):
        return True

    elif isinstance(self, eff.StateEffector):
        return False

    elif isinstance(self, eff.NotEffector):
        return _contains_boolean_contingency_effectors(self.expr)

    elif isinstance(self, eff.BinaryEffector):
        return _contains_boolean_contingency_effectors(self.left_expr) or _contains_boolean_contingency_effectors(self.right_expr)

    else:
        raise AssertionError


def _create_boolean_contingency_lookup_table(boolean_contingencies: List[ContingencyListEntry]) -> Dict[BooleanContingencyName, eff.Effector]:
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
            effector = reduce(eff.AndEffector, effector_terms)

        elif boolean_operator == BooleanOperator.op_or:
            assert len(effector_terms) > 1
            effector = reduce(eff.OrEffector, effector_terms)

        elif boolean_operator == BooleanOperator.op_not:
            assert len(effector_terms) == 1
            effector = eff.NotEffector(effector_terms[0])

        else:
            raise AssertionError

        lookup_table[current_contingency.subject] = effector

    return lookup_table


def _unary_effector_from_boolean_contingency_entry(entry: ContingencyListEntry) -> eff.Effector:
    if isinstance(entry.agent, sta.State):
        return eff.StateEffector(entry.agent)

    elif isinstance(entry.agent, BooleanContingencyName):
        return _BooleanContingencyEffector(entry.agent)

    else:
        raise AssertionError
