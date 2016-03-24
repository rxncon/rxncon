import typing as tg

import rxncon.core.state as sta
import rxncon.semantics.molecule_definition as mol
import rxncon.semantics.molecule_instance as mins
import rxncon.venntastic.sets as venn
from rxncon.semantics.molecule_definition_from_rxncon import mol_modifier_from_state_modifier


def property_set_from_mol_def_and_state_set(mol_def: mol.MoleculeDefinition, state_set: venn.Set) -> venn.Set:
    if state_set.is_equivalent_to(venn.EmptySet()):
        raise NotImplementedError

    elif state_set.is_equivalent_to(venn.UniversalSet()):
        return venn.UniversalSet()

    elif isinstance(state_set, venn.PropertySet):
        properties = _properties(mol_def, state_set.value, negate=False)
        if not properties:
            return venn.UniversalSet()
        else:
            return venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in properties], venn.Union)

    elif isinstance(state_set, venn.Complement):
        if not isinstance(state_set.expr, venn.PropertySet):
            return property_set_from_mol_def_and_state_set(mol_def, state_set.simplified_form())

        properties = _properties(mol_def, state_set.expr.value, negate=True)
        if not properties:
            return venn.UniversalSet()
        else:
            return venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in properties], venn.Union)

    elif isinstance(state_set, venn.Union):
        return venn.Union(
            property_set_from_mol_def_and_state_set(mol_def, state_set.left_expr),
            property_set_from_mol_def_and_state_set(mol_def, state_set.right_expr)
        )

    elif isinstance(state_set, venn.Intersection):
        return venn.Intersection(
            property_set_from_mol_def_and_state_set(mol_def, state_set.left_expr),
            property_set_from_mol_def_and_state_set(mol_def, state_set.right_expr)
        )


# CREATING A CONCRETE MOLECULE INSTANCE FROM THE SET OF ASSOC/MOD/LOC INSTANCES
def mol_instance_from_mol_def_and_property_set(mol_def: mol.MoleculeDefinition, property_set: venn.Set) -> mins.MoleculeInstance:
    nested_instances = property_set.to_nested_list_form()
    assert len(nested_instances) == 1

    instance_property_sets = nested_instances[0]

    instances = []
    for instance_property_set in instance_property_sets:
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
        if isinstance(instance, mins.ModificationPropertyInstance):
            mod_instances.add(instance)
        elif isinstance(instance, mins.AssociationPropertyInstance):
            assoc_instances.add(instance)
        elif isinstance(instance, mins.LocalizationPropertyInstance):
            if loc:
                raise AssertionError
            loc = instance

    return mins.MoleculeInstance(mol_def, mod_instances, assoc_instances, loc)


def mol_instance_matches_state(mol_inst: mins.MoleculeInstance, state: sta.State, negate: bool) -> bool:
    matching_properties = _properties(mol_inst.mol_def, state, negate)

    molecule_properties = []
    for x in mol_inst.association_properties:
        if x:
            molecule_properties.append(x)

    for x in mol_inst.modification_properties:
        if x:
            molecule_properties.append(x)

    if mol_inst.localization_property:
        molecule_properties.append(mol_inst.localization_property)

    return any(x in matching_properties for x in molecule_properties)


def mol_def_and_property_match_state(mol_def: mol.MoleculeDefinition, prop: mins.PropertyInstance, state: tg.Optional[sta.State], negate: bool) -> bool:

    if state is None:
        if isinstance(prop, mins.ModificationPropertyInstance):
            return prop.modifier == mins.Modifier.unmodified
        elif isinstance(prop, mins.AssociationPropertyInstance):
            return prop.occupation_status in [mins.OccupationStatus.not_occupied, mins.OccupationStatus.not_specified]
        else:
            raise NotImplementedError

    # if state is not empty and the lhs/rhs is not a part of it then it is an universal set
    if isinstance(state, sta.CovalentModificationState) and prop is not None and not mol_def.spec.is_superspecification_of(state.substrate):
        return True

    matching_properties = _properties(mol_def, state, negate)

    return prop in matching_properties


# PROTECTED HELPERS
def _properties(mol_def: mol.MoleculeDefinition, state: tg.Optional[sta.State], negate: bool) -> tg.List[mins.PropertyInstance]:

    if isinstance(state, sta.CovalentModificationState):

        matching_defs = [x for x in mol_def.modification_defs if x.spec.is_subspecification_of(state.substrate)]
        matching_instances = []
        for matching_def in matching_defs:
            if not negate:
                matching_instances.append(mins.ModificationPropertyInstance(matching_def,
                                                                            mol_modifier_from_state_modifier(state.modifier)))
            else:
                matching_instances.extend(mins.ModificationPropertyInstance(matching_def,
                                                                            mol_modifier_from_state_modifier(state.modifier)).complementary_instances())

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
                matching_instances.append(mins.AssociationPropertyInstance(matching_def,
                                                                           mol.OccupationStatus.occupied_known_partner,
                                                                           partners[0]))
            else:
                matching_instances.extend(mins.AssociationPropertyInstance(matching_def,
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
                matching_instances.append(mins.AssociationPropertyInstance(matching_def,
                                                                           mol.OccupationStatus.occupied_known_partner,
                                                                           partners[0]))
            else:
                matching_instances.extend(mins.AssociationPropertyInstance(matching_def,
                                                                           mol.OccupationStatus.occupied_known_partner,
                                                                           partners[0]).complementary_instances())

    else:
        raise NotImplementedError

    return matching_instances

