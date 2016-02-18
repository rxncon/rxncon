import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.rxncon_system as rxs
import rxncon.core.specification as spe
import rxncon.semantics.molecule_definition as mol
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.syntax.rxncon_from_string as rfs


# TESTING MOLECULE DEFINITION CREATION BY MOLECULEDEFINITIONSUPERVISOR
def test_molecule_definitions_no_contingencies():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c, b_ppi_e, b_pplus_e], [])

    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    spec_A = spe.Specification("A", None, None, None)
    spec_B = spe.Specification("B", None, None, None)
    spec_C = spe.Specification("C", None, None, None)
    spec_E = spe.Specification("E", None, None, None)

    spec_A_Bassoc = spe.Specification("A", 'Bassoc', None, None)
    spec_B_Aassoc = spe.Specification("B", 'Aassoc', None, None)
    spec_A_Cassoc = spe.Specification("A", 'Cassoc', None, None)
    spec_C_Aassoc = spe.Specification("C", 'Aassoc', None, None)

    spec_B_Eassoc = spe.Specification("B", 'Eassoc', None, None)
    spec_E_Bassoc = spe.Specification("E", 'Bassoc', None, None)
    spec_E_Bsite = spe.Specification("E", None, None, 'Bsite')

    expected_mol_def_A = mol.MoleculeDefinition(
        spec_A,
        set(),
        {mol.AssociationDefinition(spec_A_Bassoc, {spec_B_Aassoc}),
         mol.AssociationDefinition(spec_A_Cassoc, {spec_C_Aassoc})},
        None
    )

    expected_mol_def_B = mol.MoleculeDefinition(
        spec_B,
        set(),
        {mol.AssociationDefinition(spec_B_Eassoc, {spec_E_Bassoc}),
         mol.AssociationDefinition(spec_B_Aassoc, {spec_A_Bassoc})},
        None
    )

    expected_mol_def_C = mol.MoleculeDefinition(
        spec_C,
        set(),
        {mol.AssociationDefinition(spec_C_Aassoc, {spec_A_Cassoc})},
        None
    )

    expected_mol_def_E = mol.MoleculeDefinition(
        spec_E,
        {mol.ModificationDefinition(spec_E_Bsite, {mol.Modifier.unmodified, mol.Modifier.phosphorylated})},
        {mol.AssociationDefinition(spec_E_Bassoc, {spec_B_Eassoc})},
        None
    )

    assert mol_defs.molecule_definition_for_name('A') == expected_mol_def_A
    assert mol_defs.molecule_definition_for_name('B') == expected_mol_def_B
    assert mol_defs.molecule_definition_for_name('C') == expected_mol_def_C
    assert mol_defs.molecule_definition_for_name('E') == expected_mol_def_E


def test_molecule_definitions_with_contingencies():
    # Contains the same molecules (having the same definitions) as the previous tests, but this RxnCon system includes
    # contingencies.
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')

    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    b_ppi_e = rfs.reaction_from_string('B_ppi_E')
    b_dash_e = rfs.state_from_string('B--E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')

    cont_b_pplus_e = con.Contingency(b_pplus_e, con.ContingencyType.requirement, eff.StateEffector(b_dash_e))  # B_p+_E; ! B--E
    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c, b_ppi_e, b_pplus_e], [cont_b_pplus_e, cont_b_dash_e, cont_a_ppi_b])

    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    spec_A = spe.Specification("A", None, None, None)
    spec_B = spe.Specification("B", None, None, None)
    spec_C = spe.Specification("C", None, None, None)
    spec_E = spe.Specification("E", None, None, None)

    spec_A_Bassoc = spe.Specification("A", 'Bassoc', None, None)
    spec_B_Aassoc = spe.Specification("B", 'Aassoc', None, None)

    spec_A_Cassoc = spe.Specification("A", 'Cassoc', None, None)
    spec_C_Aassoc = spe.Specification("C", 'Aassoc', None, None)

    spec_B_Eassoc = spe.Specification("B", 'Eassoc', None, None)
    spec_E_Bassoc = spe.Specification("E", 'Bassoc', None, None)
    spec_E_p = spe.Specification("E", None, None, 'Bsite')

    expected_mol_def_A = mol.MoleculeDefinition(
        spec_A,
        set(),
        {mol.AssociationDefinition(spec_A_Bassoc, {spec_B_Aassoc}),
         mol.AssociationDefinition(spec_A_Cassoc, {spec_C_Aassoc})},
        None
    )

    expected_mol_def_B = mol.MoleculeDefinition(
        spec_B,
        set(),
        {mol.AssociationDefinition(spec_B_Eassoc, {spec_E_Bassoc}),
         mol.AssociationDefinition(spec_B_Aassoc, {spec_A_Bassoc})},
        None
    )

    expected_mol_def_C = mol.MoleculeDefinition(
        spec_C,
        set(),
        {mol.AssociationDefinition(spec_C_Aassoc, {spec_A_Cassoc})},
        None
    )

    expected_mol_def_E = mol.MoleculeDefinition(
        spec_E,
        {mol.ModificationDefinition(spec_E_p, {mol.Modifier.unmodified, mol.Modifier.phosphorylated})},
        {mol.AssociationDefinition(spec_E_Bassoc, {spec_B_Eassoc})},
        None
    )

    assert mol_defs.molecule_definition_for_name('A') == expected_mol_def_A
    assert mol_defs.molecule_definition_for_name('B') == expected_mol_def_B
    assert mol_defs.molecule_definition_for_name('C') == expected_mol_def_C
    assert mol_defs.molecule_definition_for_name('E') == expected_mol_def_E


def test_molecule_definitions_multiple_kinases_same_modification_at_same_residue():
    a_pplus_b = rfs.reaction_from_string('A_p+_B_[(x)]')
    c_pplus_b = rfs.reaction_from_string('C_p+_B_[(x)]')

    rxncon = rxs.RxnConSystem([a_pplus_b, c_pplus_b], [])

    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    spec_A = spe.Specification('A', None, None, None)
    spec_B = spe.Specification("B", None, None, None)
    spec_C = spe.Specification('C', None, None, None)

    spec_B_p = spe.Specification("B", None, None, 'x')

    expected_mol_def_A = mol.MoleculeDefinition(
        spec_A,
        set(),
        set(),
        None
    )

    expected_mol_def_B = mol.MoleculeDefinition(
        spec_B,
        {mol.ModificationDefinition(spec_B_p, {mol.Modifier.unmodified, mol.Modifier.phosphorylated})},
        set(),
        None
    )

    expected_mol_def_C = mol.MoleculeDefinition(
        spec_C,
        set(),
        set(),
        None
    )

    assert mol_defs.molecule_definition_for_name('A') == expected_mol_def_A
    assert mol_defs.molecule_definition_for_name('B') == expected_mol_def_B
    assert mol_defs.molecule_definition_for_name('C') == expected_mol_def_C


def test_molecule_definitions_multiple_kinases_different_modifications_at_same_residue():
    a_pplus_b = rfs.reaction_from_string('A_p+_B_[(x)]')
    c_ubplus_b = rfs.reaction_from_string('C_ub+_B_[(x)]')

    rxncon = rxs.RxnConSystem([a_pplus_b, c_ubplus_b], [])

    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    spec_A = spe.Specification('A', None, None, None)
    spec_B = spe.Specification("B", None, None, None)
    spec_C = spe.Specification('C', None, None, None)
    spec_B_p = spe.Specification("B", None, None, 'x')

    expected_mol_def_A = mol.MoleculeDefinition(
        spec_A,
        set(),
        set(),
        None
    )

    expected_mol_def_B = mol.MoleculeDefinition(
        spec_B,
        {mol.ModificationDefinition(spec_B_p, {mol.Modifier.unmodified, mol.Modifier.phosphorylated, mol.Modifier.ubiquitinated})},
        set(),
        None
    )

    expected_mol_def_C = mol.MoleculeDefinition(
        spec_C,
        set(),
        set(),
        None
    )

    assert mol_defs.molecule_definition_for_name('A') == expected_mol_def_A
    assert mol_defs.molecule_definition_for_name('B') == expected_mol_def_B
    assert mol_defs.molecule_definition_for_name('C') == expected_mol_def_C


def test_molecule_definitions_multiple_partners_binding_same_domain():
    a_pplus_b = rfs.reaction_from_string('A_ppi_B_[x]')
    c_pplus_b = rfs.reaction_from_string('C_ppi_B_[x]')

    rxncon = rxs.RxnConSystem([a_pplus_b, c_pplus_b], [])

    mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon)

    spec_A = spe.Specification('A', None, None, None)
    spec_B = spe.Specification("B", None, None, None)
    spec_C = spe.Specification('C', None, None, None)

    spec_B_x = spe.Specification("B", "x", None, None)
    spec_A_Bassoc = spe.Specification("A", "Bassoc", None, None)
    spec_C_Bassoc = spe.Specification("C", "Bassoc", None, None)

    expected_mol_def_A = mol.MoleculeDefinition(
        spec_A,
        set(),
        {mol.AssociationDefinition(spec_A_Bassoc, {spec_B_x})},
        None
    )

    expected_mol_def_B = mol.MoleculeDefinition(
        spec_B,
        set(),
        {mol.AssociationDefinition(spec_B_x, {spec_A_Bassoc, spec_C_Bassoc})},
        None
    )

    expected_mol_def_C = mol.MoleculeDefinition(
        spec_C,
        set(),
        {mol.AssociationDefinition(spec_C_Bassoc, {spec_B_x})},
        None
    )

    assert mol_defs.molecule_definition_for_name('A') == expected_mol_def_A
    assert mol_defs.molecule_definition_for_name('B') == expected_mol_def_B
    assert mol_defs.molecule_definition_for_name('C') == expected_mol_def_C

