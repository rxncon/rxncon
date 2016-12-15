import re
from typing import Dict, List, Union, Tuple, Optional
from collections import defaultdict

from typecheck import typecheck

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.effector import StateEffector, NotEffector, OrEffector, Effector, AndEffector, \
    BOOLEAN_CONTINGENCY_REGEX, BooleanOperator, BooleanContingencyName, QualSpec, qual_spec_from_str, StructEquivalences
from rxncon.core.reaction import Reaction
from rxncon.core.reaction import reaction_from_str
from rxncon.core.state import state_from_str, State


class ContingencyListEntry:
    def __init__(self, subject: Union[Reaction, BooleanContingencyName],
                 verb: Union[BooleanOperator, ContingencyType],
                 object: Union[State, BooleanContingencyName, Tuple[QualSpec, QualSpec]]):
        self.subject = subject
        self.verb    = verb
        self.object  = object

    def __eq__(self, other: 'ContingencyListEntry') -> bool:
        return self.subject == other.subject and self.verb == other.verb and self.object == other.object

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "ContingencyListEntry<{}, {}, {}>".format(self.subject, self.verb, self.object)

    @property
    def is_boolean_effector_entry(self) -> bool:
        return self.verb in (BooleanOperator(BooleanOperator.op_and), BooleanOperator(BooleanOperator.op_or),
                             BooleanOperator(BooleanOperator.op_not))

    @property
    def is_boolean_equivalence_entry(self) -> bool:
        return self.verb == BooleanOperator(BooleanOperator.op_eqv)

    @property
    def is_reaction_entry(self) -> bool:
        return isinstance(self.subject, Reaction)


class BooleanContingencyNameWithEquivs(BooleanContingencyName):
    def __init__(self, name: str, equivs: StructEquivalences):
        super().__init__(name)
        self.equivs = equivs


def contingency_list_entry_from_strs(subject_str, verb_str, object_str) -> ContingencyListEntry:
    subject_str, verb_str, object_str = subject_str.strip(), verb_str.lower().strip(), object_str.strip()

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
        equivs_dict = {int(i): qual_spec_from_str(qual_spec_str).with_prepended_namespace([BooleanContingencyName(name)])
                       for i, qual_spec_str in equivs_strs}
        equivs = StructEquivalences()
        for index, spec in enumerate(subject.components_lhs):
            try:
                equivs.add_equivalence(QualSpec([], spec.with_struct_index(index)), equivs_dict[index])
            except KeyError:
                pass

        object = BooleanContingencyNameWithEquivs(name, equivs)
    else:
        try:
            object = state_from_str(object_str)
        except SyntaxError:
            strs = [x.strip() for x in object_str.split(',')]
            object = (qual_spec_from_str(strs[0]).with_prepended_namespace([subject]),
                      qual_spec_from_str(strs[1]).with_prepended_namespace([subject]))

    return ContingencyListEntry(subject, verb, object)


@typecheck
def contingencies_from_contingency_list_entries(entries: List[ContingencyListEntry]) -> List[Contingency]:
    contingencies = []

    boolean_entries  = [x for x in entries if x.is_boolean_effector_entry]
    equiv_entries    = [x for x in entries if x.is_boolean_equivalence_entry]
    reaction_entries = [x for x in entries if x.is_reaction_entry]

    effectors        = _create_boolean_contingency_to_effector(boolean_entries)
    equivalences     = _create_boolean_contingency_to_equivalences(equiv_entries)

    while reaction_entries:
        entry = reaction_entries.pop()
        contingencies.append(Contingency(entry.subject,
                                         ContingencyType(entry.verb),
                                         _unary_effector_from_boolean_contingency_entry(entry)))

    Effector.dereference = _dereference_boolean_contingency_effectors
    Effector.contains_booleans = _contains_boolean_contingency_effectors

    while any(x.effector.contains_booleans() for x in contingencies):
        for contingency in contingencies:
            contingency.effector.dereference(effectors, equivalences)

    del Effector.dereference
    del Effector.contains_booleans

    return contingencies


class _BooleanContingencyEffector(Effector):
    def __init__(self, expr: BooleanContingencyName, equivs: Optional[StructEquivalences]=None):
        self.expr   = expr
        if not equivs:
            self.equivs = StructEquivalences()
        else:
            self.equivs = equivs

    @typecheck
    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, _BooleanContingencyEffector) and self.expr == other.expr

    def states(self):
        return [self.expr]


def _dereference_boolean_contingency_effectors(self: Effector,
                                               effector_table: Dict[BooleanContingencyName, Effector],
                                               equivalence_table: Dict[BooleanContingencyName, List[Tuple[QualSpec, QualSpec]]]):
    if isinstance(self, _BooleanContingencyEffector):
        name   = self.expr.name
        equivs = self.equivs
        self.__class__ = effector_table[self.expr].__class__
        self.__dict__  = effector_table[self.expr].__dict__
        self.name = name
        try:
            self.equivs.merge_with(equivs, [])
        except NameError:
            self.equivs = equivs
    elif isinstance(self, StateEffector):
        pass
    elif isinstance(self, NotEffector):
        _dereference_boolean_contingency_effectors(self.expr, effector_table, equivalence_table)
    elif isinstance(self, OrEffector) or isinstance(self, AndEffector):
        try:
            equivs_list = equivalence_table[BooleanContingencyName(self.name)]

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


@typecheck
def _create_boolean_contingency_to_effector(boolean_contingencies: List[ContingencyListEntry]) \
        -> Dict[BooleanContingencyName, Effector]:
    lookup_table = {}

    if not boolean_contingencies:
        return lookup_table

    assert all(x.is_boolean_effector_entry for x in boolean_contingencies)

    while boolean_contingencies:
        current_contingency = boolean_contingencies[0]
        current_contingencies = [x for x in boolean_contingencies if x.subject == current_contingency.subject]
        boolean_contingencies = [x for x in boolean_contingencies if x.subject != current_contingency.subject]

        boolean_operator = BooleanOperator(current_contingency.verb)
        assert all([BooleanOperator(x.verb) == boolean_operator for x in current_contingencies])

        effector_terms = [_unary_effector_from_boolean_contingency_entry(x) for x in current_contingencies]

        if boolean_operator == BooleanOperator.op_and:
            assert len(effector_terms) > 1
            effector = AndEffector(*effector_terms)
        elif boolean_operator == BooleanOperator.op_or:
            assert len(effector_terms) > 1
            effector = OrEffector(*effector_terms)
        elif boolean_operator == BooleanOperator.op_not:
            assert len(effector_terms) == 1
            effector = NotEffector(effector_terms[0])
        else:
            raise AssertionError

        lookup_table[current_contingency.subject] = effector

    return lookup_table


def _create_boolean_contingency_to_equivalences(equivalence_contingencies: List[ContingencyListEntry]) \
        -> Dict[BooleanContingencyName, List[Tuple[QualSpec, QualSpec]]]:
    lookup_table = defaultdict(list)

    if not equivalence_contingencies:
        return lookup_table

    assert all(x.is_boolean_equivalence_entry for x in equivalence_contingencies)

    for contingency in equivalence_contingencies:
        lookup_table[contingency.subject].append(contingency.object)

    return lookup_table


@typecheck
def _unary_effector_from_boolean_contingency_entry(entry: ContingencyListEntry) -> Effector:
    if isinstance(entry.object, State):
        return StateEffector(entry.object)
    elif isinstance(entry.object, BooleanContingencyNameWithEquivs):
        return _BooleanContingencyEffector(entry.object, entry.object.equivs)
    elif isinstance(entry.object, BooleanContingencyName):
        return _BooleanContingencyEffector(entry.object)
    else:
        raise AssertionError
