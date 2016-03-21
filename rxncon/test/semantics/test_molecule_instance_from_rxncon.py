
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.venntastic.sets as venn
import rxncon.core.rxncon_system as rxs
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.specification as spe
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.semantics.molecule_instance_from_rxncon as mir
import rxncon.semantics.molecule_definition as mdf
import rxncon.semantics.molecule_instance as mins
import rxncon.simulation.rule_based.rbm_from_rxncon as rfr


def test_set_of_instances_from_molecule_def_and_set_of_states_for_ppi_no_contingency():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    rxncon = rxs.RxnConSystem([a_ppi_b], [])
    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    strict_cont_state_set = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_b))
    assert strict_cont_state_set.is_equivalent_to(venn.UniversalSet())

    # Now test the molecule instances that get created by taking into account the strict contingencies. Note that this
    # does not take into account the source contingency (i.e. A not bound to B).
    strict_instances_set = mir.property_set_from_mol_def_and_state_set(
        mol_defs.mol_def_for_name('A'),
        strict_cont_state_set
    )
    assert strict_instances_set.is_equivalent_to(venn.UniversalSet())

    strict_instances_set = mir.property_set_from_mol_def_and_state_set(
        mol_defs.mol_def_for_name('B'),
        strict_cont_state_set
    )
    assert strict_instances_set.is_equivalent_to(venn.UniversalSet())


def test_set_of_instances_from_molecule_def_and_set_of_states_for_ppi_and_requirement_contingency():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    mol_def_A = mol_defs.mol_def_for_name('A')

    strict_cont_state_set = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_b))
    strict_instances_set = mir.property_set_from_mol_def_and_state_set(mol_def_A, strict_cont_state_set)

    assoc_defs_A_to_C = [assoc_def for assoc_def in mol_def_A.association_defs if assoc_def.spec.domain == "Cassoc"]
    assert len(assoc_defs_A_to_C) == 1
    assoc_def = assoc_defs_A_to_C[0]

    assert isinstance(assoc_def, mdf.AssociationPropertyDefinition)
    assert assoc_def.spec == spe.Specification('A', 'Cassoc', None, None)

    expected_assoc_instance = mins.AssociationPropertyInstance(assoc_def,
                                                               mins.OccupationStatus.occupied_known_partner,
                                                               spe.Specification('C', 'Aassoc', None, None))

    assert strict_instances_set.is_equivalent_to(venn.PropertySet(expected_assoc_instance))


def test_set_of_instances_from_molecule_def_and_set_of_states_for_ppi_and_inhibition_contingency():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    mol_def_A = mol_defs.mol_def_for_name('A')

    strict_cont_state_set = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_b))
    strict_instances_set = mir.property_set_from_mol_def_and_state_set(mol_def_A, strict_cont_state_set)

    assoc_defs_A_to_C = [assoc_def for assoc_def in mol_def_A.association_defs if assoc_def.spec.domain == "Cassoc"]
    assert len(assoc_defs_A_to_C) == 1
    assoc_def = assoc_defs_A_to_C[0]

    assert isinstance(assoc_def, mdf.AssociationPropertyDefinition)
    assert assoc_def.spec == spe.Specification('A', 'Cassoc', None, None)

    expected_assoc_instance = mins.AssociationPropertyInstance(assoc_def,
                                                               mdf.OccupationStatus.not_occupied,
                                                               None)

    assert strict_instances_set.is_equivalent_to(venn.PropertySet(expected_assoc_instance))


def test_set_of_instances_from_complex_system():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')
    b_ppi_e = rfs.reaction_from_string('B_ppi_E')
    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string("E-{P}")

    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(e_pplus))  # B_ppi_E; ! E-{P}

    rxncon = rxs.RxnConSystem([b_ppi_e, a_ppi_b, b_pplus_e], [cont_e_pplus, cont_b_dash_e])
    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    # TEST MOLECULE A
    actual_A_set_of_instances = mir.property_set_from_mol_def_and_state_set(
        mol_defs.mol_def_for_name('A'),
        rfr.state_set_from_contingencies([cont_b_dash_e, cont_e_pplus])
    )

    assoc_def_A_to_B = mdf.AssociationPropertyDefinition(spe.Specification('A', 'Bassoc', None, None),
                                                         {spe.Specification('B', 'Aassoc', None, None)})
    expected_A_set_of_instances = venn.PropertySet(
        mins.AssociationPropertyInstance(assoc_def_A_to_B,
                                         mdf.OccupationStatus.occupied_known_partner,
                                         spe.Specification('B', 'Aassoc', None, None)))

    assert actual_A_set_of_instances.is_equivalent_to(expected_A_set_of_instances)

    # TEST MOLECULE B
    actual_B_set_of_instances = mir.property_set_from_mol_def_and_state_set(
        mol_defs.mol_def_for_name('B'),
        rfr.state_set_from_contingencies([cont_b_dash_e, cont_e_pplus])
    )

    assoc_def_B_to_A = mdf.AssociationPropertyDefinition(spe.Specification('B', 'Aassoc', None, None),
                                                         {spe.Specification('A', 'Bassoc', None, None)})
    expected_B_set_of_instances = venn.PropertySet(
        mins.AssociationPropertyInstance(assoc_def_B_to_A,
                                         mdf.OccupationStatus.occupied_known_partner,
                                         spe.Specification('A', 'Bassoc', None, None)))

    assert actual_B_set_of_instances.is_equivalent_to(expected_B_set_of_instances)

    # TEST MOLECULE E
    actual_E_set_of_instances = mir.property_set_from_mol_def_and_state_set(
        mol_defs.mol_def_for_name('E'),
        rfr.state_set_from_contingencies([cont_b_dash_e, cont_e_pplus])
    )

    mod_def_E = mdf.ModificationPropertyDefinition(
        spe.Specification('E', None, None, 'Bsite'),
        {mdf.Modifier.unmodified, mdf.Modifier.phosphorylated}
    )

    expected_E_set_of_instances = venn.PropertySet(mins.ModificationPropertyInstance(mod_def_E,
                                                                                     mdf.Modifier.phosphorylated))

    assert actual_E_set_of_instances.is_equivalent_to(expected_E_set_of_instances)


def test_molecule_instance_matches_state():
    spec_A = spe.Specification("A", None, None, None)
    spec_B = spe.Specification("B", None, None, None)

    spec_A_Bassoc = spe.Specification("A", 'Bassoc', None, None)
    spec_B_Aassoc = spe.Specification("B", 'Aassoc', None, None)

    spec_A_psite = spe.Specification('A', None, None, 'psite')

    mol_def_A = mdf.MoleculeDefinition(
        spec_A,
        {mdf.ModificationPropertyDefinition(spec_A_psite, {mdf.Modifier.unmodified, mdf.Modifier.phosphorylated})},
        {mdf.AssociationPropertyDefinition(spec_A_Bassoc, {spec_B_Aassoc})},
        None
    )

    assoc_inst_A_bound = mins.AssociationPropertyInstance(mdf.AssociationPropertyDefinition(spec_A_Bassoc, {spec_B_Aassoc}),
                                                          mdf.OccupationStatus.occupied_known_partner,
                                                          spec_B_Aassoc)

    assoc_inst_A_free = mins.AssociationPropertyInstance(mdf.AssociationPropertyDefinition(spec_A_Bassoc, {spec_B_Aassoc}),
                                                         mdf.OccupationStatus.not_occupied,
                                                         None)

    mod_inst_A_phos = mins.ModificationPropertyInstance(mdf.ModificationPropertyDefinition(spec_A_psite, {mdf.Modifier.unmodified, mdf.Modifier.phosphorylated}),
                                                        mdf.Modifier.phosphorylated)

    mod_inst_A_unphos = mins.ModificationPropertyInstance(mdf.ModificationPropertyDefinition(spec_A_psite, {mdf.Modifier.unmodified, mdf.Modifier.phosphorylated}),
                                                          mdf.Modifier.unmodified)

    # A BOUND, UNPHOSPHORYLATED
    mol_inst_A_bound_unphos = mins.MoleculeInstance(mol_def_A, {mod_inst_A_unphos}, {assoc_inst_A_bound}, None)
    assert mir.mol_instance_matches_state(mol_inst_A_bound_unphos,
                                          rfs.state_from_string('A--B'),
                                          False)

    assert not mir.mol_instance_matches_state(mol_inst_A_bound_unphos,
                                              rfs.state_from_string('A--B'),
                                              True)

    assert not mir.mol_instance_matches_state(mol_inst_A_bound_unphos,
                                              rfs.state_from_string('A-{p}'),
                                              False)

    assert mir.mol_instance_matches_state(mol_inst_A_bound_unphos,
                                          rfs.state_from_string('A-{p}'),
                                          True)

    # A FREE, UNPHOSPHORYLATED
    mol_inst_A_unbound_unphos = mins.MoleculeInstance(mol_def_A, {mod_inst_A_unphos}, {assoc_inst_A_free}, None)
    assert not mir.mol_instance_matches_state(mol_inst_A_unbound_unphos,
                                              rfs.state_from_string('A--B'),
                                              False)

    assert mir.mol_instance_matches_state(mol_inst_A_unbound_unphos,
                                          rfs.state_from_string('A--B'),
                                          True)

    assert not mir.mol_instance_matches_state(mol_inst_A_unbound_unphos,
                                              rfs.state_from_string('A-{p}'),
                                              False)

    assert mir.mol_instance_matches_state(mol_inst_A_unbound_unphos,
                                          rfs.state_from_string('A-{p}'),
                                          True)

    # A BOUND, PHOSPHORYLATED
    mol_inst_A_bound_phos = mins.MoleculeInstance(mol_def_A, {mod_inst_A_phos}, {assoc_inst_A_bound}, None)
    assert mir.mol_instance_matches_state(mol_inst_A_bound_phos,
                                          rfs.state_from_string('A--B'),
                                          False)

    assert not mir.mol_instance_matches_state(mol_inst_A_bound_phos,
                                              rfs.state_from_string('A--B'),
                                              True)

    assert mir.mol_instance_matches_state(mol_inst_A_bound_phos,
                                          rfs.state_from_string('A-{p}'),
                                          False)

    assert not mir.mol_instance_matches_state(mol_inst_A_bound_phos,
                                              rfs.state_from_string('A-{p}'),
                                              True)

    # A BOUND, PHOSPHORYLATED
    mol_inst_A_free_phos = mins.MoleculeInstance(mol_def_A, {mod_inst_A_phos}, {assoc_inst_A_free}, None)
    assert not mir.mol_instance_matches_state(mol_inst_A_free_phos,
                                              rfs.state_from_string('A--B'),
                                              False)

    assert mir.mol_instance_matches_state(mol_inst_A_free_phos,
                                          rfs.state_from_string('A--B'),
                                          True)

    assert mir.mol_instance_matches_state(mol_inst_A_free_phos,
                                          rfs.state_from_string('A-{p}'),
                                          False)

    assert not mir.mol_instance_matches_state(mol_inst_A_free_phos,
                                              rfs.state_from_string('A-{p}'),
                                              True)



