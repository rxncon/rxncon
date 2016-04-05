from typing import Dict, List, Tuple, Optional
from itertools import product
from collections import defaultdict
from rxncon.util.utils import compose

from rxncon.venntastic.sets import Set, Complement, Union, Intersection, PropertySet, nested_expression_from_list_and_binary_op
from rxncon.core.specification import Specification
from rxncon.core.reaction import Reaction, ReactionClass, Verb
from rxncon.core.state import State, CovalentModificationState, InterProteinInteractionState, IntraProteinInteractionState, \
    TranslocationState, InputState, SynthesisDegradationState, StateModifier
from rxncon.semantics.molecule_definition import MoleculeDefinition, Modifier
from rxncon.semantics.molecule_definition_from_rxncon import mod_domain_spec_from_state_and_reaction, ass_domain_specs_from_state
from rxncon.semantics.molecule_instance import MoleculeInstance, ModificationPropertyInstance, AssociationPropertyInstance, \
    OccupationStatus


STATE_TO_MOLECULE_INSTANCE_FUNCTIONS = {
    CovalentModificationState:    '_molecule_instance_set_from_mod_state',
    InterProteinInteractionState: '_molecule_instance_set_from_ppi_state',
    IntraProteinInteractionState: '_molecule_instance_set_from_ipi_state',
    TranslocationState:           '_molecule_instance_set_from_loc_state',
    InputState:                   '_molecule_instance_set_from_inp_state'
}


REACTION_TO_MOLECULE_INSTANCE_PAIRS_FUNCTIONS = {
    Verb.phosphorylation:             '_molecule_instance_set_pair_from_pplus_reaction',
    Verb.protein_protein_interaction: '_molecule_instance_set_pair_from_ppi_reaction'
}


def mol_instance_set_from_state_set(mol_defs: Dict[Specification, MoleculeDefinition], state_set: Set) -> Set:
    def _mol_instance_set_with_complements_from_state_set(mol_defs: Dict[Specification, MoleculeDefinition], state_set: Set) -> Set:
        if isinstance(state_set, PropertySet):
            assert isinstance(state_set.value, State)
            return globals()[STATE_TO_MOLECULE_INSTANCE_FUNCTIONS[type(state_set.value)]](mol_defs, state_set.value)
        elif isinstance(state_set, Complement):
            return Complement(_mol_instance_set_with_complements_from_state_set(mol_defs, state_set.expr))
        elif isinstance(state_set, Intersection):
            return Intersection(_mol_instance_set_with_complements_from_state_set(mol_defs, state_set.left_expr),
                                _mol_instance_set_with_complements_from_state_set(mol_defs, state_set.right_expr))
        elif isinstance(state_set, Union):
            return Union(_mol_instance_set_with_complements_from_state_set(mol_defs, state_set.left_expr),
                         _mol_instance_set_with_complements_from_state_set(mol_defs, state_set.right_expr))
        else:
            raise NotImplemented

    def _implode_mol_instance_set(mol_instance_set: Set) -> Set:
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
                assert isinstance(mol.value, MoleculeInstance)
                instances[mol.value.mol_def] += [mol.value]

            imploded_mols = [_imploded_mol_instances(v) for k, v in instances.items() if _imploded_mol_instances(v)]
            cleaned_terms.append(nested_expression_from_list_and_binary_op([PropertySet(x) for x in imploded_mols], Intersection))

        return nested_expression_from_list_and_binary_op(cleaned_terms, Union)

    def _expanded_complements(mol_instance_set: Set) -> Set:
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
                    raise AssertionError

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
            assert isinstance(mol_instance_set.value, MoleculeInstance)
            return mol_instance_set
        else:
            raise NotImplemented

    mols_from_states = lambda x: _mol_instance_set_with_complements_from_state_set(mol_defs, x)

    full_mapping = compose(
        _implode_mol_instance_set,
        _expanded_complements,
        mols_from_states
    )

    return full_mapping(state_set).simplified_form()


def mol_instance_set_pair_from_reaction(mol_defs: Dict[Specification, MoleculeDefinition], reaction: Reaction) -> Tuple[Set, Set]:
    subject_mol_def = mol_defs[Specification(reaction.subject.name, None, None, None)]
    object_mol_def = mol_defs[Specification(reaction.object.name, None, None, None)]

    return globals()[REACTION_TO_MOLECULE_INSTANCE_PAIRS_FUNCTIONS[reaction.verb]](reaction, subject_mol_def, object_mol_def)


### MAPPING STATES TO MOLECULE INSTANCES ###
def _molecule_instance_set_from_mod_state(mol_defs: Dict[Specification, MoleculeDefinition], state: CovalentModificationState) -> Set:
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


    mol_def = mol_defs[Specification(state.substrate.name, None, None, None)]

    mod_defs = [mod_def for mod_def in mol_def.modification_defs if
                mod_def.spec.is_subspecification_of(state.substrate)]
    modifier = _mol_modifier_from_state_modifier(state.modifier)
    mol_instances = [
        PropertySet(MoleculeInstance(mol_def, {ModificationPropertyInstance(x, modifier)}, set(), None))
        for x in mod_defs]

    return nested_expression_from_list_and_binary_op(mol_instances, Union)

def _molecule_instance_set_from_ppi_state(mol_defs: Dict[Specification, MoleculeDefinition], state: InterProteinInteractionState) -> Set:
    first_mol_def = mol_defs[Specification(state.first_component.name, None, None, None)]
    second_mol_def = mol_defs[Specification(state.second_component.name, None, None, None)]

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

def _molecule_instance_set_from_ipi_state(mol_defs: Dict[Specification, MoleculeDefinition], state: IntraProteinInteractionState) -> Set:
    pass

def _molecule_instance_set_from_loc_state(mol_defs: Dict[Specification, MoleculeDefinition], state: TranslocationState) -> Set:
    pass

def _molecule_instance_set_from_inp_state(mol_defs: Dict[Specification, MoleculeDefinition], state: InputState) -> Set:
    pass

### MAPPING REACTIONS TO MOLECULE INSTANCE PAIRS ###
def _molecule_instance_set_pair_from_pplus_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[Set, Set]:
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

def _molecule_instance_set_pair_from_ppi_reaction(reaction: Reaction, subject_mol_def: MoleculeDefinition, object_mol_def: MoleculeDefinition) -> Tuple[Set, Set]:
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

