from collections import defaultdict
from copy import deepcopy
from enum import Enum
from typing import List, Optional, Dict, Tuple
from typecheck import typecheck

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.contingency import ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.spec import Spec, LocusResolution, BondSpec
from rxncon.core.state import StateDef, State
from rxncon.core.reaction import Reaction, ReactionTerm, ReactionTerm, BondReactionTerm
from rxncon.util.utils import all_eq
from rxncon.venntastic.sets import Set, ValueSet, Union, Intersection, Complement, MultiIntersection, MultiUnion


class RuleBasedModel:
    @typecheck
    def __init__(self, molecule_defs: List[MolDef], rules: List['Rule'],
                 parameters: List['Parameter'], initial_conditions: List['InitialCondition']):
        self.mol_defs = molecule_defs
        self.rules = rules
        self.parameters = parameters
        self.initial_conditions = initial_conditions

        self._validate()

    def set_parameter_value(self, parameter_name, parameter_value):
        # @todo
        pass

    def set_initial_condition(self, molecule, value):
        # @todo
        pass

    def _validate(self):
        for initial_condition in self.initial_conditions:
            if initial_condition.molecule.mol_def not in self.mol_defs:
                raise ValueError('Initial condition {0} refers to unknown molecule def {1}.'
                                 .format(initial_condition, initial_condition.molecule.mol_def))

        for rule in self.rules:
            for molecule in rule.molecules:
                if molecule not in self.mol_defs:
                    raise ValueError('Rule {0} contains molecule def {1}, which is absent in the model'
                                     .format(rule, molecule))


class Parameter:
    @typecheck
    def __init__(self, name: str, value: Optional[str]):
        self.name = name
        self.value = value

    @typecheck
    def __eq__(self, other: 'Parameter') -> bool:
        return self.name == other.name and self.value == other.value

    def __hash__(self):
        return hash(self.name)

    def __lt__(self, other: 'Parameter'):
        if self.name != other.name:
            return self.name < other.name
        if not self.value:
            return True
        elif not other.value:
            return False
        else:
            return self.value < other.value

    def __repr__(self):
        return str(self)

    def __str__(self) -> bool:
        return 'Parameter: {0} = {1}'.format(self.name, self.value)


class InitialCondition:
    def __init__(self, molecule: Mol, value):
        self.molecule = molecule
        self.value = value

    def __eq__(self, other: 'InitialCondition'):
        return self.molecule == other.molecule and self.value == other.value

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'InitialCondition: {0} = {1}'.format(self.molecule, self.value)


class MutualExclusivityError(Exception):
    pass


class InconsistentStructureError(Exception):
    pass


class InsufficientStructureError(Exception):
    pass


class MolDef:
    @typecheck
    def __init__(self, component: Spec, valid_states: Dict[Tuple[Spec, StateDef], List[State]]):
        assert component.has_resolution(LocusResolution.component)
        assert not component.struct_index

        self.component    = component
        self.valid_states = valid_states
        self._valid_states_nested = defaultdict(dict)
        for (spec, state_def), states in self.valid_states.items():
            self._valid_states_nested[spec][state_def] = states

        self._validate()

    def __str__(self) -> str:
        return 'MolDef<{0}>, valid_states: <{1}>'.format(str(self.component), self.valid_states)

    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def valid_states_by_spec(self, spec: Spec) -> Dict[StateDef, List[State]]:
        return self._valid_states_nested[spec]

    @typecheck
    def complementary_states(self, spec: Spec, state: State) -> List[State]:
        return [x for x in self.valid_states[(spec, state.definition)] if x != state]

    def _validate(self):
        for (spec, state_def), states in self.valid_states.items():
            assert all(spec in state.specs for state in states)
            assert all_eq([state_def] + [state.definition for state in states])
            assert not spec.struct_index
            assert not any(state.is_structured for state in states)


class Mol:
    @typecheck
    def __init__(self, mol_def: MolDef, states: Dict[Tuple[Spec, StateDef], State]):
        self.mol_def   = mol_def
        self.states    = states
        self.component = deepcopy(mol_def.component)
        self._mutation_post_process()

    def __str__(self):
        return 'Mol<{0}>, states: <{1}>'.format(str(self.mol_def.component), str(self.states))

    def __repr__(self):
        return str(self)

    def add_state(self, mol_spec: Spec, state: State):
        assert not mol_spec.struct_index
        assert isinstance(state.target, Spec)

        if (mol_spec, state.definition) in self.states.keys() and self.states[(mol_spec, state.definition)] != state:
            raise MutualExclusivityError

        self.states[(mol_spec, state.definition)] = state
        self._mutation_post_process()

    def _mutation_post_process(self):
        self._determine_struct_index()
        self._validate()

    def _determine_struct_index(self):
        UNIVERSAL_SET = {'all'}

        possible_indices = UNIVERSAL_SET

        for state in self.states.values():
            indices = {component.struct_index for component in state.components if component.is_non_struct_equiv_to(self.component)
                       and component.struct_index is not None}

            if indices and possible_indices == UNIVERSAL_SET:
                possible_indices = set(indices)
            elif indices:
                possible_indices = possible_indices & indices

        if possible_indices == set():
            raise InconsistentStructureError('Inconsistent structure annotation on Mol {}'.format(str(self)))
        elif len(possible_indices) > 1:
            raise InsufficientStructureError('Insufficient structure annotation on Mol {}'.format(str(self)))
        elif possible_indices != UNIVERSAL_SET and len(possible_indices) == 1:
            self.component.struct_index = possible_indices.pop()

    def _validate(self):
        for (spec, state_def), state in self.states.items():
            assert state.to_non_structured_state() in self.mol_def.valid_states[(spec, state_def)]


@typecheck
def mol_def_for_component(rxncon_sys: RxnConSystem, component: Spec) -> MolDef:
    return MolDef(component, rxncon_sys.states_for_component_grouped(component))


class Complex:
    def __init__(self, mols: List[Mol], bonds: List[BondSpec]):
        self.mols  = {mol.component: mol for mol in mols}
        self.bonds = bonds

    def add_mol(self, mol: Mol):
        if self.contains_component(mol.component):
            raise AssertionError

        self.mols[mol.component] = mol

    def add_bond(self, bond: BondSpec):
        self.bonds.append(bond)

    def contains_component(self, component: Spec) -> bool:
        return component in self.mols.keys()

    @property
    def is_closed(self):
        return True


class Arrow(Enum):
    single_headed = 0
    double_headed = 1


class Rule:
    def __init__(self, lhs: List[Complex], rhs: List[Complex], arrow: Arrow, rates: List[Parameter]):
        self.lhs, self.rhs, self.arrow, self.rates = lhs, rhs, arrow, rates


STATE_TO_COMPLEX = {

}





class RuleBuilder:
    def __init__(self, rxncon_sys: RxnConSystem):
        self.rxncon_sys = rxncon_sys
        self.mol_defs   = {}
        self._determine_mol_defs()

    def rules_from_reaction(self, reaction: Reaction) -> List[Rule]:
        rules = []

        lhs_center = [self._complex_from_rxn_term(term) for term in reaction.terms_lhs]
        rhs_center = [self._complex_from_rxn_term(term) for term in reaction.terms_rhs]

        overlapping_contexts = self._contexts_for_reaction(reaction)
        finished_contexts = []

        while overlapping_contexts:
            # Each of these 'overlapping' contexts originates straight from the contingencies. Therefore these
            # must always be connected to the reactants. The problem with non-connectedness arises when we make them disjunct.
            current_context   = overlapping_contexts.pop()
            disjunct_contexts = self._make_disjunct(current_context, finished_contexts)

            while disjunct_contexts:
                disjunct_context = disjunct_contexts.pop()
                if not self._connected(disjunct_context, lhs_center):
                    # If we fail, try again later.
                    disjunct_contexts.insert(0, disjunct_context)
                else:
                    finished_contexts.append(disjunct_context)
                    rules.append(Rule(self._make_complexes(lhs_center, disjunct_context),
                                      self._make_complexes(rhs_center, disjunct_context),
                                      self._arrow_for_reaction(reaction),
                                      self._rates_for_reaction(reaction)))

        return rules

    def _connected(self, context: Set, center: List[Complex]) -> bool:
        bonds = [bond for bond in context.values

    @typecheck
    def _make_disjunct(self, context: Set, finished_contexts: List[Set]) -> List[Set]:
        disjunct_contexts = []

        for s in Intersection(context, Complement(MultiUnion(*finished_contexts))).to_intersection_terms():
            no_complements     = self._expand_complements(s)
            disjunct_contexts += no_complements.to_intersection_terms()

        return disjunct_contexts

    def _determine_mol_defs(self):
        for component in self.rxncon_sys.components:
            self.mol_defs[component] = mol_def_for_component(self.rxncon_sys, component)

    @typecheck
    def _expand_complements(self, s: Set) -> Set:
        if isinstance(s, ValueSet):
            return s
        elif isinstance(s, Intersection):
            return Intersection(self._expand_complements(s.left_expr), self._expand_complements(s.right_expr))
        elif isinstance(s, Union):
            return Union(self._expand_complements(s.left_expr), self._expand_complements(s.right_expr))
        elif isinstance(s, Complement):
            if isinstance(s.expr, Complement):
                return self._expand_complements(s.expr.expr)
            elif isinstance(s.expr, Union):
                return Intersection(self._expand_complements(Complement(s.expr.left_expr)),
                                    self._expand_complements(Complement(s.expr.right_expr)))
            elif isinstance(s.expr, Intersection):
                return Union(self._expand_complements(Complement(s.expr.left_expr)),
                             self._expand_complements(Complement(s.expr.right_expr)))
            elif isinstance(s.expr, ValueSet):
                assert isinstance(s.expr.value, State)
                terms = []
                for spec in s.expr.value.specs:
                    component   = spec.to_non_struct_spec().to_component_spec()
                    complements = self.mol_defs[component].complementary_states(spec, s.expr.value)
                    terms.append(MultiUnion(*[ValueSet(x) for x in complements]))

                return MultiIntersection(*terms)

    @typecheck
    def _complex_from_rxn_term(self, rxn_term: ReactionTerm) -> Complex:
        if isinstance(rxn_term, ReactionTerm):
            mol = Mol(self.mol_defs[rxn_term.spec], {})
            for state in rxn_term.states:
                mol.add_state(state.target, state)

            return Complex([mol], rxn_term.bonds)
        # @todo JCR We do not handle states on bonds yet.

    @typecheck
    def _contexts_for_reaction(self, reaction: Reaction) -> List[List[State]]:
        contingencies = self.rxncon_sys.s_contingencies_for_reaction(reaction)

        effectors = [x.effector for x in contingencies if x.type == ContingencyType.requirement] + \
                    [NotEffector(x.effector) for x in contingencies if x.type == ContingencyType.inhibition]

        context = MultiIntersection(*[_venn_set_from_effector(x) for x in effectors])
        context = self._expand_complements(context)

        return context.to_intersection_terms()

def _venn_set_from_effector(effector: Effector) -> Set:
    if isinstance(effector, StateEffector):
        return ValueSet(effector.expr)
    elif isinstance(effector, AndEffector):
        return Intersection(
            _venn_set_from_effector(effector.left_expr),
            _venn_set_from_effector(effector.right_expr)
        )
    elif isinstance(effector, OrEffector):
        return Union(
            _venn_set_from_effector(effector.left_expr),
            _venn_set_from_effector(effector.right_expr)
        )
    elif isinstance(effector, NotEffector):
        return Complement(_venn_set_from_effector(effector.expr))
    else:
        raise Exception('Unknown Effector {}'.format(effector))

