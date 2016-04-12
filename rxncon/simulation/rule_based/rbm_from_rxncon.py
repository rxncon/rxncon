from typing import List, Set, Tuple
import itertools as itt


from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction, Directionality as RxnDirectionality, Influence as RxnInfluence
from rxncon.core.state import State
from rxncon.core.contingency import ContingencyType
from rxncon.core.effector import Effector, AndEffector, NotEffector, OrEffector, StateEffector
from rxncon.semantics.molecule_instance import MoleculeInstance
from rxncon.simulation.rule_based.rule_based_model import Reactant, ComplexReactant, MoleculeReactant
from rxncon.simulation.rule_based.molecule_from_rxncon import imploded_mol_instance_set
from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, Rule, Arrow, Parameter
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys
from rxncon.simulation.rule_based.molecule_from_rxncon import mol_instance_set_from_state_set, mol_instance_set_pair_from_reaction
from rxncon.venntastic.sets import Set as VennSet, PropertySet, Complement, Union, Intersection, nested_expression_from_list_and_binary_op


def rbm_from_rxncon_sys(rxconsys: RxnConSystem) -> RuleBasedModel:
    rules = set()

    for reaction in rxconsys.reactions:
        rules = rules.union(rules_from_reaction(rxconsys, reaction))

    mol_defs = mol_defs_from_rxncon_sys(rxconsys)

    return RuleBasedModel(set(mol_defs.values()), rules, set(), set())


def rules_from_reaction(rxconsys: RxnConSystem, reaction: Reaction) -> Set[Rule]:
    def get_arrow(reaction):
        if reaction.influence == RxnInfluence.bidirectional:
            return Arrow.reversible
        elif reaction.influence in [RxnInfluence.positive, RxnInfluence.negative, RxnInfluence.transfer]:
            return Arrow.irreversible
        else:
            raise AssertionError

    def get_rates(reaction, qcc: _QuantitativeContingencyConfiguration):
        rate_suffix = '_{}'.format(str(reaction))
        if str(qcc):
            rate_suffix += '_{}'.format(str(qcc))

        if reaction.influence in [RxnInfluence.positive, RxnInfluence.negative, RxnInfluence.transfer]:
            return {
                Parameter('k' + rate_suffix, None)
            }
        elif reaction.influence == RxnInfluence.bidirectional:
            return {
                Parameter('kf' + rate_suffix, None),
                Parameter('kr' + rate_suffix, None)
            }

    quant_cont_configs = _quant_contingency_configs_from_reaction(rxconsys, reaction)
    mol_defs = mol_defs_from_rxncon_sys(rxconsys)

    rules = set()

    for qcc in quant_cont_configs:
        background_molecule_sets = \
            mol_instance_set_from_state_set(
                mol_defs,
                Intersection(qcc.to_state_set(), _strict_contingency_state_set_from_reaction(rxconsys, reaction))
            ).to_union_list_form()

        lhs_reacting_molecule_set, rhs_reacting_molecule_set = mol_instance_set_pair_from_reaction(mol_defs, reaction)

        for background_molecule_set in background_molecule_sets:
            lhs_reactants = _reactants_from_molecule_sets(lhs_reacting_molecule_set, background_molecule_set)
            rhs_reactants = _reactants_from_molecule_sets(rhs_reacting_molecule_set, background_molecule_set)

            rules.add(Rule(lhs_reactants, rhs_reactants, get_arrow(reaction), get_rates(reaction, qcc)))

    return rules


def _reactants_from_molecule_sets(reacting_molecule_set: VennSet, background_molecule_set: VennSet) -> Set[Reactant]:
    def get_molecules():
        molecules = imploded_mol_instance_set(Intersection(reacting_molecule_set, background_molecule_set)).to_nested_list_form()

        assert len(molecules) == 1
        assert all(isinstance(x, PropertySet) for x in molecules[0])
        molecules = [x.value for x in molecules[0]]

        reacting_molecules = [molecule for molecule in molecules if molecule.mol_def in
                              (x.value.mol_def for x in reacting_molecule_set.to_nested_list_form()[0])]

        background_molecules = [molecule for molecule in molecules if molecule not in reacting_molecules]
        return sorted(reacting_molecules), sorted(background_molecules)

    reacting_molecules, background_molecules = get_molecules()

    reactants = []
    complexes_in_formation = []

    for molecule in reacting_molecules:
        if not molecule.bindings:
            reactants.append(MoleculeReactant(molecule))
            continue

        still_single = True
        for complex in complexes_in_formation:
            if complex.can_bind_molecule_instance(molecule):
                complex.add_molecule_instance(molecule)
                still_single = False
                break

        if still_single:
            new_complex = _ComplexInFormation()
            new_complex.add_molecule_instance(molecule)
            complexes_in_formation.append(new_complex)

    failed_binding_attempts = 0
    while background_molecules:
        molecule = background_molecules.pop(0)
        if not molecule.bindings:
            # Spurious solution, background molecule not connected.
            return set()

        still_single = True
        for complex in complexes_in_formation:
            if complex.can_bind_molecule_instance(molecule):
                complex.add_molecule_instance(molecule)
                still_single = False
                failed_binding_attempts = 0
                break

        if still_single:
            background_molecules.append(molecule)
            failed_binding_attempts += 1
            if failed_binding_attempts > len(background_molecules):
                return set()

    reactants += [cif.to_complex_reactant() for cif in complexes_in_formation]

    return set(reactants)


class _ComplexInFormation:
    def __init__(self):
        self.molecules = set()
        self.realized_bindings = set()
        self.potential_bindings = set()

    def can_bind_molecule_instance(self, mol: MoleculeInstance):
        return not self.molecules or any(mol_binding in self.potential_bindings for mol_binding in mol.bindings)

    def add_molecule_instance(self, mol: MoleculeInstance):
        assert self.can_bind_molecule_instance(mol)

        if not self.molecules:
            self.molecules.add(mol)
            self.potential_bindings = mol.bindings
        else:
            # Greedily fill the bindings.
            self.molecules.add(mol)
            for new_binding in sorted(mol.bindings):
                if new_binding in self.potential_bindings:
                    self.potential_bindings.remove(new_binding)
                    self.realized_bindings.add(new_binding)
                else:
                    self.potential_bindings.add(new_binding)

    def to_complex_reactant(self):
        assert not self.potential_bindings
        return ComplexReactant(self.molecules, self.realized_bindings)


def _strict_contingency_state_set_from_reaction(rxnconsys: RxnConSystem, reaction: Reaction) -> VennSet:
    def _state_set_from_effector(effector: Effector) -> VennSet:
        if isinstance(effector, StateEffector):
            return PropertySet(effector.expr)
        if isinstance(effector, NotEffector):
            return Complement(_state_set_from_effector(effector.expr))
        elif isinstance(effector, AndEffector):
            return Intersection(_state_set_from_effector(effector.left_expr),
                                _state_set_from_effector(effector.right_expr))
        elif isinstance(effector, OrEffector):
            return Union(_state_set_from_effector(effector.left_expr), _state_set_from_effector(effector.right_expr))
        else:
            raise AssertionError

    required = [_state_set_from_effector(cont.effector) for cont in rxnconsys.strict_contingencies_for_reaction(reaction)
                if cont.type == ContingencyType.requirement]
    verboten = [_state_set_from_effector(cont.effector) for cont in rxnconsys.strict_contingencies_for_reaction(reaction)
                if cont.type == ContingencyType.inhibition]

    return Intersection(nested_expression_from_list_and_binary_op(required, Intersection),
                        Complement(nested_expression_from_list_and_binary_op(verboten, Union)))


def _quant_contingency_configs_from_reaction(rxnconsys: RxnConSystem, reaction: Reaction) -> Set['_QuantitativeContingencyConfiguration']:
    def _true_false_combinations(xs: Set) -> List[Tuple[Set, Set]]:
        result = []
        for n in range(len(xs)):
            true_combis = itt.combinations(xs, n)
            for trues in true_combis:
                result.append((set(trues), {f for f in xs if f not in trues}))

        return result

    states = set()
    for contingency in rxnconsys.quantitative_contingencies_for_reaction(reaction):
        states = states.union(set(contingency.effector.states))

    combis = _true_false_combinations(states)

    if combis:
        return {_QuantitativeContingencyConfiguration(combi[0], combi[1]) for combi in combis}
    else:
        return {_QuantitativeContingencyConfiguration(set(), set())}


class _QuantitativeContingencyConfiguration:
    def __init__(self, present_states: Set[State], absent_states: Set[State]):
        self.present_states = present_states
        self.absent_states = absent_states

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return '!{0}'.join(sorted(self.present_states)) + 'x{0}'.join(sorted(self.absent_states))

    def to_state_set(self) -> VennSet:
        return Intersection(
            nested_expression_from_list_and_binary_op([PropertySet(x) for x in self.present_states], Intersection),
            nested_expression_from_list_and_binary_op([Complement(PropertySet(x)) for x in self.absent_states], Intersection)
        )
