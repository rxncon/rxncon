import typing as tg

import rxncon.core.state as sta
import rxncon.semantics.molecule_definition as mol
import rxncon.semantics.molecule_instance as mins
import rxncon.venntastic.sets as venn
from rxncon.semantics.molecule_definition_from_rxncon import molecule_modifier_from_state_modifier


def set_of_instances_from_molecule_def_and_set_of_states(mol_def: mol.MoleculeDefinition, set_of_states: venn.Set) -> venn.Set:
    if set_of_states.is_equivalent_to(venn.EmptySet()):
        raise NotImplementedError

    elif set_of_states.is_equivalent_to(venn.UniversalSet()):
        return venn.UniversalSet()

    elif isinstance(set_of_states, venn.PropertySet):
        instances = _instances(mol_def, set_of_states.value, negate=False)
        if not instances:
            return venn.UniversalSet()
        else:
            return venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in instances], venn.Union)

    elif isinstance(set_of_states, venn.Complement):
        if not isinstance(set_of_states.expr, venn.PropertySet):
            return set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.simplified_form())

        instances = _instances(mol_def, set_of_states.expr.value, negate=True)
        if not instances:
            return venn.UniversalSet()
        else:
            return venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in instances], venn.Union)

    elif isinstance(set_of_states, venn.Union):
        return venn.Union(
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.left_expr),
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.right_expr)
        )

    elif isinstance(set_of_states, venn.Intersection):
        return venn.Intersection(
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.left_expr),
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.right_expr)
        )


# CREATING A CONCRETE MOLECULE INSTANCE FROM THE SET OF ASSOC/MOD/LOC INSTANCES
def molecule_instance_from_molecule_def_and_set_of_instances(mol_def: mol.MoleculeDefinition, set_of_instances: venn.Set) -> mins.MoleculeInstance:
    nested_instances = set_of_instances.to_nested_list_form()
    assert len(nested_instances) == 1

    instance_property_sets = nested_instances[0]

    instances = []
    for instance_property_set in instance_property_sets:
        print(instance_property_set)
        if isinstance(instance_property_set, venn.PropertySet):
            instances.append(instance_property_set.value)
        elif isinstance(instance_property_set, venn.EmptySet):
            break
        else:
            raise AssertionError

    assoc_instances = set()
    mod_instances = set()
    loc = None

    for instance in instances:
        if isinstance(instance, rxncon.semantics.molecule_instance.ModificationInstance):
            mod_instances.add(instance)
        elif isinstance(instance, rxncon.semantics.molecule_instance.AssociationInstance):
            assoc_instances.add(instance)
        elif isinstance(instance, rxncon.semantics.molecule_instance.LocalizationInstance):
            if loc:
                raise AssertionError
            loc = instance

    return rxncon.semantics.molecule_instance.MoleculeInstance(mol_def, mod_instances, assoc_instances, loc)


def molecule_instance_matches_state(mol_inst: mins.MoleculeInstance, state: sta.State, negate: bool) -> bool:
    matching_instances = _instances(mol_inst.molecule_def, state, negate)

    molecule_instances = []
    for x in mol_inst.association_instances:
        if x:
            molecule_instances.append(x)

    for x in mol_inst.modification_instances:
        if x:
            molecule_instances.append(x)

    if mol_inst.localization_instance:
        molecule_instances.append(mol_inst.localization_instance)

    return any(x in matching_instances for x in molecule_instances)


# PROTECTED HELPERS
def _instances(mol_def: mol.MoleculeDefinition, state: sta.State, negate: bool) -> tg.List[mins.Instance]:
    if isinstance(state, sta.CovalentModificationState):
        matching_defs = [x for x in mol_def.modification_defs if x.spec.is_subspecification_of(state.substrate)]
        matching_instances = []
        for matching_def in matching_defs:
            if not negate:
                matching_instances.append(mins.ModificationInstance(matching_def,
                                                                    molecule_modifier_from_state_modifier(state.modifier)))
            else:
                matching_instances.extend(mins.ModificationInstance(matching_def,
                                                                    molecule_modifier_from_state_modifier(state.modifier)).complementary_instances())

        return matching_instances

    elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
        first_defs = [assoc_def for assoc_def in mol_def.association_defs
                      for valid_partner_spec in assoc_def.valid_partners
                      if assoc_def.spec.is_subspecification_of(state.first_component) and
                      valid_partner_spec.is_subspecification_of(state.second_component)]
        assert len(first_defs) <= 1
        matching_instances = []
        for matching_def in first_defs:
            partners = [spec for spec in matching_def.valid_partners if spec.is_subspecification_of(state.second_component)]
            assert len(partners) == 1

            if not negate:
                matching_instances.append(mins.AssociationInstance(matching_def,
                                                                   mol.OccupationStatus.occupied_known_partner,
                                                                   partners[0]))
            else:
                matching_instances.extend(mins.AssociationInstance(matching_def,
                                                                   mol.OccupationStatus.occupied_known_partner,
                                                                   partners[0]).complementary_instances())

        second_defs = [assoc_def for assoc_def in mol_def.association_defs
                       for valid_partner_spec in assoc_def.valid_partners
                       if assoc_def.spec.is_subspecification_of(state.second_component) and
                       valid_partner_spec.is_subspecification_of(state.first_component)]
        assert len(second_defs) <= 1
        for matching_def in second_defs:
            partners = [spec for spec in matching_def.valid_partners if spec.is_subspecification_of(state.first_component)]
            assert len(partners) == 1

            if not negate:
                matching_instances.append(mins.AssociationInstance(matching_def,
                                                                   mol.OccupationStatus.occupied_known_partner,
                                                                   partners[0]))
            else:
                matching_instances.extend(mins.AssociationInstance(matching_def,
                                                                   mol.OccupationStatus.occupied_known_partner,
                                                                   partners[0]).complementary_instances())

    else:
        raise NotImplementedError

    return matching_instances

