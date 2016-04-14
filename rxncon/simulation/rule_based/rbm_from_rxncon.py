from collections import defaultdict
from itertools import product
from typing import List, Set, Tuple, Dict, Optional
import itertools as itt


from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction, Directionality as RxnDirectionality, Influence as RxnInfluence, Verb
from rxncon.core.specification import Specification
from rxncon.core.state import State, CovalentModificationState, InterProteinInteractionState, \
    IntraProteinInteractionState, TranslocationState, InputState, StateModifier
from rxncon.core.contingency import ContingencyType
from rxncon.core.effector import Effector, AndEffector, NotEffector, OrEffector, StateEffector
from rxncon.semantics.molecule_definition import MoleculeDefinition, Modifier, OccupationStatus
from rxncon.semantics.molecule_instance import MoleculeInstance, ModificationPropertyInstance, \
    AssociationPropertyInstance
from rxncon.simulation.rule_based.rule_based_model import Reactant, ComplexReactant, MoleculeReactant
from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, Rule, Arrow, Parameter
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys, \
    mod_domain_spec_from_state_and_reaction, ass_domain_specs_from_state
from rxncon.util.utils import compose
from rxncon.venntastic.sets import Set as VennSet, PropertySet, Complement, Union, Intersection, nested_expression_from_list_and_binary_op, \
    UniversalSet, EmptySet, gram_schmidt_disjunctify


STATE_TO_MOLECULE_INSTANCE_FUNCTIONS = {
    CovalentModificationState:    '_molecule_instance_set_from_mod_state',
    InterProteinInteractionState: '_molecule_instance_set_from_ppi_state',
    IntraProteinInteractionState: '_molecule_instance_set_from_ipi_state',
    TranslocationState:           '_molecule_instance_set_from_loc_state',
    InputState:                   '_molecule_instance_set_from_inp_state'
}
REACTION_TO_MOLECULE_INSTANCE_PAIRS_FUNCTIONS = {
    Verb.phosphorylation:             '_molecule_instance_set_pair_from_pplus_reaction',
    Verb.dephosphorylation:           '_molecule_instance_set_pair_from_pminus_reaction',
    Verb.phosphotransfer:             '_molecule_instance_set_pair_from_pt_reaction',
    Verb.ubiquitination:              '_molecule_instance_set_pair_from_ubplus_reaction',
    Verb.deubiquitination:            '_molecule_instance_set_pair_from_ubminus_reaction',
    Verb.guanine_nucleotide_exchange: '_molecule_instance_set_pair_from_gef_reaction',
    Verb.gtpase_activation:           '_molecule_instance_set_pair_from_gap_reaction',
    Verb.proteolytic_cleavage:        '_molecule_instance_set_pair_from_cut_reaction',
    Verb.protein_protein_interaction: '_molecule_instance_set_pair_from_ppi_reaction',
}



def rbm_from_rxncon_sys(rxconsys: RxnConSystem) -> RuleBasedModel:
    rules = set()

    for reaction in rxconsys.reactions:
        rules = rules.union(rules_from_reaction(rxconsys, reaction))

    mol_defs = mol_defs_from_rxncon_sys(rxconsys)

    return RuleBasedModel(set(mol_defs.values()), rules, set(), set())


def rules_from_reaction(rxconsys: RxnConSystem, reaction: Reaction) -> Set[Rule]:
    def get_arrow(reaction):
        if reaction.directionality == RxnDirectionality.bidirectional:
            return Arrow.reversible
        elif reaction.directionality == RxnDirectionality.unidirectional:
            return Arrow.irreversible
        else:
            raise AssertionError

    def get_rates(reaction, qcc: _QuantitativeContingencyConfiguration):
        rate_suffix = '_{}'.format(str(reaction))
        if str(qcc):
            rate_suffix += '_{}'.format(str(qcc))

        if reaction.directionality == RxnDirectionality.unidirectional:
            return {
                Parameter('k' + rate_suffix, None)
            }
        elif reaction.directionality == RxnDirectionality.bidirectional:
            return {
                Parameter('kf' + rate_suffix, None),
                Parameter('kr' + rate_suffix, None)
            }

    quant_cont_configs = _quant_contingency_configs_from_reaction(rxconsys, reaction)
    mol_defs = mol_defs_from_rxncon_sys(rxconsys)
    lhs_reacting_molecule_set, rhs_reacting_molecule_set = mol_instance_set_pair_from_reaction(mol_defs, reaction)

    rules = set()

    for qcc in quant_cont_configs:
        background_molecule_sets = \
            mol_instance_set_from_state_set(
                mol_defs,
                Intersection(qcc.to_state_set(), _strict_contingency_state_set_from_reaction(rxconsys, reaction)),
                lhs_reacting_molecule_set
            ).to_union_list_form()

        for background_molecule_set in background_molecule_sets:
            lhs_reactants = reactants_from_molecule_sets(lhs_reacting_molecule_set, background_molecule_set)
            rhs_reactants = reactants_from_molecule_sets(rhs_reacting_molecule_set, background_molecule_set)

            if lhs_reactants and rhs_reactants:
                rules.add(Rule(lhs_reactants, rhs_reactants, get_arrow(reaction), get_rates(reaction, qcc)))

    return rules


def reactants_from_molecule_sets(reacting_molecule_set: VennSet, background_molecule_set: VennSet) -> Set[Reactant]:
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

    if any(not x.is_finished for x in complexes_in_formation):
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

    @property
    def is_finished(self) -> bool:
        return not self.potential_bindings

    def to_complex_reactant(self):
        assert self.is_finished
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
        for n in range(len(xs) + 1):
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
        return ''.join('!{}'.format(str(x)) for x in sorted(self.present_states)) + \
               ''.join('x{}'.format(str(x)) for x in sorted(self.absent_states))

    def to_state_set(self) -> VennSet:
        return Intersection(
            nested_expression_from_list_and_binary_op([PropertySet(x) for x in self.present_states], Intersection),
            nested_expression_from_list_and_binary_op([Complement(PropertySet(x)) for x in self.absent_states], Intersection)
        )


def mol_instance_set_from_state_set(mol_defs: Dict[Specification, MoleculeDefinition], state_set: VennSet, reacting_set: VennSet) -> VennSet:
    def _mol_instance_set_with_complements_from_state_set(mol_defs: Dict[Specification, MoleculeDefinition], state_set: VennSet) -> VennSet:
        if isinstance(state_set, PropertySet):
            if not state_set.value:
                return UniversalSet()
            assert isinstance(state_set.value, State)
            return globals()[STATE_TO_MOLECULE_INSTANCE_FUNCTIONS[type(state_set.value)]](mol_defs, state_set.value)
        elif isinstance(state_set, EmptySet):
            return EmptySet()
        elif isinstance(state_set, Complement):
            return Complement(_mol_instance_set_with_complements_from_state_set(mol_defs, state_set.expr))
        elif isinstance(state_set, Intersection):
            return Intersection(_mol_instance_set_with_complements_from_state_set(mol_defs, state_set.left_expr),
                                _mol_instance_set_with_complements_from_state_set(mol_defs, state_set.right_expr))
        elif isinstance(state_set, Union):
            return Union(_mol_instance_set_with_complements_from_state_set(mol_defs, state_set.left_expr),
                         _mol_instance_set_with_complements_from_state_set(mol_defs, state_set.right_expr))
        else:
            raise NotImplementedError

    def _expanded_complements(mol_instance_set: VennSet) -> VennSet:
        def _exploded_mol_instance(mol_instance: MoleculeInstance) -> List[MoleculeInstance]:
            exploded = [MoleculeInstance(mol_instance.mol_def, {x}, set(), None) for x in
                        mol_instance.modification_properties]
            exploded += [MoleculeInstance(mol_instance.mol_def, set(), {x}, None) for x in
                         mol_instance.association_properties]
            if mol_instance.localization_property:
                exploded += [MoleculeInstance(mol_instance.mol_def, set(), set(), mol_instance.localization_property)]

            return exploded

        def _expanded_complementary_molecule_instance_set(mol_instance_set: Complement) -> Set:
            assert isinstance(mol_instance_set.expr, PropertySet)
            assert isinstance(mol_instance_set.expr.value, MoleculeInstance)
            mol_def = mol_instance_set.expr.value.mol_def
            single_property_instances = _exploded_mol_instance(mol_instance_set.expr.value)
            expanded_instances = []

            for instance in single_property_instances:
                if instance.modification_properties:
                    assert len(instance.modification_properties) == 1
                    assert not instance.association_properties
                    assert not instance.localization_property
                    prop = list(instance.modification_properties)[0]
                    expanded_instances += [MoleculeInstance(mol_def, {x}, set(), None) for x in
                                           prop.complementary_instances()]
                elif instance.association_properties:
                    assert len(instance.association_properties) == 1
                    assert not instance.modification_properties
                    assert not instance.localization_property
                    prop = list(instance.association_properties)[0]
                    expanded_instances += [MoleculeInstance(mol_def, set(), {x}, None) for x in
                                           prop.complementary_instances()]
                elif instance.localization_property:
                    assert not instance.modification_properties
                    assert not instance.association_properties
                    expanded_instances += [MoleculeInstance(mol_def, set(), set(), instance.localization_property)]
                else:
                    raise AssertionError('The exploded MoleculeInstance should have one and exactly one PropertyInstance')

            return nested_expression_from_list_and_binary_op([PropertySet(x) for x in expanded_instances], Union)

        mol_instance_set = mol_instance_set.simplified_form()

        if isinstance(mol_instance_set, Complement):
            assert isinstance(mol_instance_set.expr, PropertySet)
            assert isinstance(mol_instance_set.expr.value, MoleculeInstance)
            return _expanded_complementary_molecule_instance_set(mol_instance_set)
        elif isinstance(mol_instance_set, Intersection):
            return Intersection(_expanded_complements(mol_instance_set.left_expr),
                                _expanded_complements(mol_instance_set.right_expr))
        elif isinstance(mol_instance_set, Union):
            return Union(_expanded_complements(mol_instance_set.left_expr),
                         _expanded_complements(mol_instance_set.right_expr))
        elif isinstance(mol_instance_set, PropertySet):
            assert not mol_instance_set.value or isinstance(mol_instance_set.value, MoleculeInstance)
            return mol_instance_set
        else:
            raise NotImplemented

    def _filter_disconnected_components(mol_instance_set: VennSet, reacting_set: VennSet) -> VennSet:
        all_solutions = mol_instance_set.to_union_list_form()

        connected_solutions = [x for x in all_solutions if reactants_from_molecule_sets(reacting_set, x)]

        return nested_expression_from_list_and_binary_op(connected_solutions, Union)

    def _disjunctified(mol_instance_set: VennSet) -> VennSet:
        return nested_expression_from_list_and_binary_op(gram_schmidt_disjunctify(mol_instance_set.to_union_list_form()), Union)

    def _sorted(mol_instance_set: VennSet) -> VennSet:
        sorted_nested_list = []
        for sublist in mol_instance_set.to_nested_list_form():
            sorted_nested_list.append(sorted(sublist))

        sorted_nested_list.sort()

        sorted_sets = []
        for sublist in sorted_nested_list:
            sorted_sets.append(nested_expression_from_list_and_binary_op(sublist, Intersection))

        return nested_expression_from_list_and_binary_op(sorted_sets, Union)

    mols_from_states = lambda x: _mol_instance_set_with_complements_from_state_set(mol_defs, x)
    filter_disconnected = lambda x: _filter_disconnected_components(x, reacting_set)

    full_mapping = compose(
        _sorted,
        imploded_mol_instance_set,
        _expanded_complements,
        _disjunctified,
        _sorted,
        filter_disconnected,
        _expanded_complements,
        mols_from_states
    )

    return full_mapping(state_set).simplified_form()


def mol_instance_set_pair_from_reaction(mol_defs: Dict[Specification, MoleculeDefinition], reaction: Reaction) -> Tuple[VennSet, VennSet]:
    subject_mol_def = mol_defs[reaction.subject.to_component_specification()]
    object_mol_def = mol_defs[reaction.object.to_component_specification()]

    try:
        return globals()[REACTION_TO_MOLECULE_INSTANCE_PAIRS_FUNCTIONS[reaction.verb]](reaction, subject_mol_def, object_mol_def)
    except KeyError:
        raise AssertionError('Reaction verb {} missing in REACTION_TO_MOLECULE_INSTANCE_PAIRS_FUNCTIONS table.'
                             .format(str(reaction.verb)))


def imploded_mol_instance_set(mol_instance_set: VennSet) -> VennSet:
    def _imploded_mol_instances(mol_instances: List[MoleculeInstance]) -> Optional[MoleculeInstance]:
        # Returns None if the List of MoleculeInstance contains internal inconsistencies, e.g.
        # if there are multiple localization properties, or if a single residue is both phosphorylated and
        # unmodified.
        mol_def = mol_instances[0].mol_def
        assert all(x.mol_def == mol_def for x in mol_instances)
        asss = set()
        mods = set()
        loc = None

        for mol_instance in mol_instances:
            [asss.add(x) for x in mol_instance.association_properties]
            [mods.add(x) for x in mol_instance.modification_properties]
            if mol_instance.localization_property:
                if loc is None:
                    loc = mol_instance.localization_property
                else:
                    return None

        if len([ass.association_def for ass in asss]) != len(set([ass.association_def for ass in asss])):
            return None

        if len([mod.modification_def for mod in mods]) != len(set([mod.modification_def for mod in mods])):
            return None

        return MoleculeInstance(mol_def, mods, asss, loc)

    nested_form = mol_instance_set.to_nested_list_form()
    cleaned_terms = []

    for term in nested_form:
        instances = defaultdict(list)
        for mol in term:
            assert isinstance(mol, PropertySet)
            assert not mol.value or isinstance(mol.value, MoleculeInstance)
            if mol.value:
                instances[mol.value.mol_def] += [mol.value]

        imploded_mols = [_imploded_mol_instances(v) for k, v in instances.items() if _imploded_mol_instances(v)]
        cleaned_terms.append(
            nested_expression_from_list_and_binary_op([PropertySet(x) for x in imploded_mols], Intersection))

    return nested_expression_from_list_and_binary_op(cleaned_terms, Union)


def _molecule_instance_set_from_mod_state(mol_defs: Dict[Specification, MoleculeDefinition], state: CovalentModificationState) -> VennSet:
    def _mol_modifier_from_state_modifier(state_modifier: StateModifier) -> Modifier:
        if state_modifier == StateModifier.unmodified:
            return Modifier.unmodified
        elif state_modifier == StateModifier.phosphor:
            return Modifier.phosphorylated
        elif state_modifier == StateModifier.ubiquitin:
            return Modifier.ubiquitinated
        elif state_modifier == StateModifier.truncated:
            return Modifier.truncated
        else:
            raise NotImplemented


    mol_def = mol_defs[state.substrate.to_component_specification()]

    mod_defs = [mod_def for mod_def in mol_def.modification_defs if
                mod_def.spec.is_subspecification_of(state.substrate)]
    modifier = _mol_modifier_from_state_modifier(state.modifier)
    mol_instances = [
        PropertySet(MoleculeInstance(mol_def, {ModificationPropertyInstance(x, modifier)}, set(), None))
        for x in mod_defs]

    return nested_expression_from_list_and_binary_op(mol_instances, Union)


def _molecule_instance_set_from_ppi_state(mol_defs: Dict[Specification, MoleculeDefinition], state: InterProteinInteractionState) -> VennSet:
    first_mol_def = mol_defs[state.first_component.to_component_specification()]
    second_mol_def = mol_defs[state.second_component.to_component_specification()]

    first_ass_defs = [ass_def for ass_def in first_mol_def.association_defs
                      if any(state.second_component.is_superspecification_of(x) for x in ass_def.valid_partners)]

    second_ass_defs = [ass_def for ass_def in second_mol_def.association_defs
                       if any(state.first_component.is_superspecification_of(x) for x in ass_def.valid_partners)]

    first_ass_instances = [AssociationPropertyInstance(ass_def, OccupationStatus.occupied_known_partner, partner)
                           for ass_def in first_ass_defs for partner in ass_def.valid_partners
                           if state.second_component.is_superspecification_of(partner)]

    second_ass_instances = [AssociationPropertyInstance(ass_def, OccupationStatus.occupied_known_partner, partner)
                            for ass_def in second_ass_defs for partner in ass_def.valid_partners
                            if state.first_component.is_superspecification_of(partner)]

    sets = [Intersection(PropertySet(MoleculeInstance(first_mol_def, set(), {pair[0]}, None)),
                         PropertySet(MoleculeInstance(second_mol_def, set(), {pair[1]}, None))) for pair in
            product(first_ass_instances, second_ass_instances)
            if pair[0].association_def.spec == pair[1].partner and pair[0].partner == pair[1].association_def.spec]

    return nested_expression_from_list_and_binary_op(sets, Union)


def _molecule_instance_set_from_ipi_state(mol_defs: Dict[Specification, MoleculeDefinition], state: IntraProteinInteractionState) -> VennSet:
    # @todo
    pass


def _molecule_instance_set_from_loc_state(mol_defs: Dict[Specification, MoleculeDefinition], state: TranslocationState) -> VennSet:
    # @todo
    pass


def _molecule_instance_set_from_inp_state(mol_defs: Dict[Specification, MoleculeDefinition], state: InputState) -> VennSet:
    # @todo
    pass


def _molecule_instance_set_pair_from_pplus_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.product, reaction)]
    assert len(mod_defs) == 1
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    phosphorylated = ModificationPropertyInstance(mod_def, Modifier.phosphorylated)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {phosphorylated}, set(), None))))


def _molecule_instance_set_pair_from_pminus_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.source, reaction)]
    assert len(mod_defs) == 1, 'Could not find phosphorylation domain for molecule {0}, for reaction {1}'.format(str(object_mol_def), str(reaction))
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    phosphorylated = ModificationPropertyInstance(mod_def, Modifier.phosphorylated)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {phosphorylated}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))))


def _molecule_instance_set_pair_from_pt_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    object_mod_defs = [x for x in object_mol_def.modification_defs
                       if x.spec == mod_domain_spec_from_state_and_reaction(reaction.product, reaction)]

    assert len(object_mod_defs) == 1, 'Could not find object phosphorylation domain for reaction {}'.format(str(reaction))
    object_mod_def = object_mod_defs[0]

    object_unmodified = ModificationPropertyInstance(object_mod_def, Modifier.unmodified)
    object_phosphorylated = ModificationPropertyInstance(object_mod_def, Modifier.phosphorylated)


    subject_mod_defs = [x for x in subject_mol_def.modification_defs
                        if x.spec == mod_domain_spec_from_state_and_reaction(reaction.source, reaction)]

    assert len(subject_mod_defs) == 1, 'Could not find subject phosphorylation domain for reaction {}'.format(str(reaction))
    subject_mod_def = subject_mod_defs[0]

    subject_unmodified = ModificationPropertyInstance(subject_mod_def, Modifier.unmodified)
    subject_phosphorylated = ModificationPropertyInstance(subject_mod_def, Modifier.phosphorylated)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, {subject_phosphorylated}, set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {object_unmodified}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, {subject_unmodified}, set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {object_phosphorylated}, set(), None))))


def _molecule_instance_set_pair_from_ubplus_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.product, reaction)]
    assert len(mod_defs) == 1
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    ubiquitinated = ModificationPropertyInstance(mod_def, Modifier.ubiquitinated)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {ubiquitinated}, set(), None))))


def _molecule_instance_set_pair_from_ubminus_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.source, reaction)]
    assert len(mod_defs) == 1, 'Could not find ubiquination domain for molecule {0}, for reaction {1}'.format(str(object_mol_def), str(reaction))
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    ubiquitinated = ModificationPropertyInstance(mod_def, Modifier.ubiquitinated)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {ubiquitinated}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))))


def _molecule_instance_set_pair_from_gef_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.product, reaction)]
    assert len(mod_defs) == 1
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    gtp = ModificationPropertyInstance(mod_def, Modifier.guanosintriphosphat)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {gtp}, set(), None))))


def _molecule_instance_set_pair_from_gap_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.source, reaction)]
    assert len(mod_defs) == 1, 'Could not find gtp domain for molecule {0}, for reaction {1}'.format(str(object_mol_def), str(reaction))
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    gtp = ModificationPropertyInstance(mod_def, Modifier.guanosintriphosphat)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {gtp}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))))


def _molecule_instance_set_pair_from_cut_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    mod_defs = [x for x in object_mol_def.modification_defs
                if x.spec == mod_domain_spec_from_state_and_reaction(reaction.product, reaction)]
    assert len(mod_defs) == 1
    mod_def = mod_defs[0]

    unmodified = ModificationPropertyInstance(mod_def, Modifier.unmodified)
    cut = ModificationPropertyInstance(mod_def, Modifier.truncated)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {unmodified}, set(), None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), set(), None)),
                         PropertySet(MoleculeInstance(object_mol_def, {cut}, set(), None))))


def _molecule_instance_set_pair_from_ppi_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[VennSet, VennSet]:
    first_ass_spec, second_ass_spec = ass_domain_specs_from_state(reaction.product)

    first_ass_defs  = [x for x in subject_mol_def.association_defs if x.spec == first_ass_spec and second_ass_spec in x.valid_partners]
    second_ass_defs = [x for x in object_mol_def.association_defs if x.spec == second_ass_spec and first_ass_spec in x.valid_partners]
    assert len(first_ass_defs) == 1
    assert len(second_ass_defs) == 1

    first_ass_def, second_ass_def = first_ass_defs[0], second_ass_defs[0]
    first_free = AssociationPropertyInstance(first_ass_def, OccupationStatus.not_occupied, None)
    first_bound = AssociationPropertyInstance(first_ass_def, OccupationStatus.occupied_known_partner, second_ass_def.spec)

    second_free = AssociationPropertyInstance(second_ass_def, OccupationStatus.not_occupied, None)
    second_bound = AssociationPropertyInstance(second_ass_def, OccupationStatus.occupied_known_partner, first_ass_def.spec)

    return (Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), {first_free}, None)),
                         PropertySet(MoleculeInstance(object_mol_def, set(), {second_free}, None))),
            Intersection(PropertySet(MoleculeInstance(subject_mol_def, set(), {first_bound}, None)),
                         PropertySet(MoleculeInstance(object_mol_def, set(), {second_bound}, None))))