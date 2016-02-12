import rxncon.semantics.molecule_from_rxncon as mfr
import rxncon.semantics.molecule as mol
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.rxncon_system as rxs
import rxncon.core.specification as spe
import rxncon.core.contingency as con
import rxncon.core.effector as eff

def test_MoleculeDefinitionSupervisor():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c, b_ppi_e, b_pplus_e], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    universal_specA = spe.Specification("A", None, None, None)
    universal_specB = spe.Specification("B", None, None, None)
    universal_specC = spe.Specification("C", None, None, None)
    universal_specE = spe.Specification("E", None, None, None)

    specificationAiB = spe.Specification("A", 'Bassoc', None, None)
    specificationBiA = spe.Specification("B", 'Aassoc', None, None)

    specificationAiC = spe.Specification("A", 'Cassoc', None, None)
    specificationCiA = spe.Specification("C", 'Aassoc', None, None)

    specificationBiE = spe.Specification("B", 'Eassoc', None, None)
    specificationEiB = spe.Specification("E", 'Bassoc', None, None)
    specificationEp = spe.Specification("E", None, None, 'Bsite')


    expected_mol_def_A = mol.MoleculeDefinition(universal_specA, set(), {mol.AssociationDefinition(specificationAiB, {specificationBiA}),
                                                             mol.AssociationDefinition(specificationAiC, {specificationCiA})},
                                                mol.LocalizationDefinition(set()))
    expected_mol_def_B = mol.MoleculeDefinition(universal_specB, set(), {mol.AssociationDefinition(specificationBiE, {specificationEiB}),
                                                                     mol.AssociationDefinition(specificationBiA, {specificationAiB})},
                                                mol.LocalizationDefinition(set()))
    expected_mol_def_C = mol.MoleculeDefinition(universal_specC, set(), {mol.AssociationDefinition(specificationCiA, {specificationAiC})},
                                                mol.LocalizationDefinition(set()))
    expected_mol_def_E = mol.MoleculeDefinition(universal_specE, {mol.ModificationDefinition(specificationEp,
                                                                                            {mol.Modifier.unmodified, mol.Modifier.phosphorylated})},
                                                                 {mol.AssociationDefinition(specificationEiB, {specificationBiE})},
                                                                 mol.LocalizationDefinition(set()))


    assert len(mol_defs.molecule_definition_for_name("A").association_defs) == 2
    assert not mol_defs.molecule_definition_for_name("A").modification_defs
    assert mol_defs.molecule_definition_for_name("A").localization_def == expected_mol_def_A.localization_def
    assert list(mol_defs.molecule_definition_for_name("A").association_defs)[0] != list(mol_defs.molecule_definition_for_name("A").association_defs)[1]
    assert list(mol_defs.molecule_definition_for_name("A").association_defs)[0] in [list(expected_mol_def_A.association_defs)[0], list(expected_mol_def_A.association_defs)[1]]
    assert list(mol_defs.molecule_definition_for_name("A").association_defs)[1] in [list(expected_mol_def_A.association_defs)[0], list(expected_mol_def_A.association_defs)[1]]

    assert not mol_defs.molecule_definition_for_name("B").modification_defs
    assert len(mol_defs.molecule_definition_for_name("B").association_defs) == 2
    assert mol_defs.molecule_definition_for_name("B").localization_def == expected_mol_def_B.localization_def
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[0] != list(mol_defs.molecule_definition_for_name("B").association_defs)[1]
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[0] in [list(expected_mol_def_B.association_defs)[0], list(expected_mol_def_B.association_defs)[1]]
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[1] in [list(expected_mol_def_B.association_defs)[0], list(expected_mol_def_B.association_defs)[1]]

    assert mol_defs.molecule_definition_for_name("C") == expected_mol_def_C

    assert mol_defs.molecule_definition_for_name("E").modification_defs == expected_mol_def_E.modification_defs
    assert mol_defs.molecule_definition_for_name("E").localization_def == expected_mol_def_E.localization_def
    assert len(mol_defs.molecule_definition_for_name("E").association_defs) == 1
    assert list(mol_defs.molecule_definition_for_name("E").association_defs)[0] == list(expected_mol_def_E.association_defs)[0]


def test_molecule_defs_from_rxncon_with_contingencies():
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

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    universal_specA = spe.Specification("A", None, None, None)
    universal_specB = spe.Specification("B", None, None, None)
    universal_specC = spe.Specification("C", None, None, None)
    universal_specE = spe.Specification("E", None, None, None)

    specificationAiB = spe.Specification("A", 'Bassoc', None, None)
    specificationBiA = spe.Specification("B", 'Aassoc', None, None)

    specificationAiC = spe.Specification("A", 'Cassoc', None, None)
    specificationCiA = spe.Specification("C", 'Aassoc', None, None)

    specificationBiE = spe.Specification("B", 'Eassoc', None, None)
    specificationEiB = spe.Specification("E", 'Bassoc', None, None)
    specificationEp = spe.Specification("E", None, None, 'Bsite')


    expected_mol_def_A = mol.MoleculeDefinition(universal_specA, set(), {mol.AssociationDefinition(specificationAiB, {specificationBiA}),
                                                             mol.AssociationDefinition(specificationAiC, {specificationCiA})},
                                                mol.LocalizationDefinition(set()))
    expected_mol_def_B = mol.MoleculeDefinition(universal_specB, set(), {mol.AssociationDefinition(specificationBiE, {specificationEiB}),
                                                                     mol.AssociationDefinition(specificationBiA, {specificationAiB})},
                                                mol.LocalizationDefinition(set()))
    expected_mol_def_C = mol.MoleculeDefinition(universal_specC, set(), {mol.AssociationDefinition(specificationCiA, {specificationAiC})},
                                                mol.LocalizationDefinition(set()))
    expected_mol_def_E = mol.MoleculeDefinition(universal_specE, {mol.ModificationDefinition(specificationEp,
                                                                                            {mol.Modifier.unmodified, mol.Modifier.phosphorylated})},
                                                                 {mol.AssociationDefinition(specificationEiB, {specificationBiE})},
                                                                 mol.LocalizationDefinition(set()))

    assert len(mol_defs.molecule_definition_for_name("A").association_defs) == 2
    assert not mol_defs.molecule_definition_for_name("A").modification_defs
    assert mol_defs.molecule_definition_for_name("A").localization_def == expected_mol_def_A.localization_def
    assert list(mol_defs.molecule_definition_for_name("A").association_defs)[0] != list(mol_defs.molecule_definition_for_name("A").association_defs)[1]
    assert list(mol_defs.molecule_definition_for_name("A").association_defs)[0] in [list(expected_mol_def_A.association_defs)[0], list(expected_mol_def_A.association_defs)[1]]
    assert list(mol_defs.molecule_definition_for_name("A").association_defs)[1] in [list(expected_mol_def_A.association_defs)[0], list(expected_mol_def_A.association_defs)[1]]

    assert not mol_defs.molecule_definition_for_name("B").modification_defs
    assert len(mol_defs.molecule_definition_for_name("B").association_defs) == 2
    assert mol_defs.molecule_definition_for_name("B").localization_def == expected_mol_def_B.localization_def
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[0] != list(mol_defs.molecule_definition_for_name("B").association_defs)[1]
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[0] in [list(expected_mol_def_B.association_defs)[0], list(expected_mol_def_B.association_defs)[1]]
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[1] in [list(expected_mol_def_B.association_defs)[0], list(expected_mol_def_B.association_defs)[1]]

    assert mol_defs.molecule_definition_for_name("C") == expected_mol_def_C

    assert mol_defs.molecule_definition_for_name("E").modification_defs == expected_mol_def_E.modification_defs
    assert mol_defs.molecule_definition_for_name("E").localization_def == expected_mol_def_E.localization_def
    assert len(mol_defs.molecule_definition_for_name("E").association_defs) == 1
    assert list(mol_defs.molecule_definition_for_name("E").association_defs)[0] == list(expected_mol_def_E.association_defs)[0]

def test_molecule_defs_from_rxncon_modifiation_at_same_residue():
    a_pplus_b = rfs.reaction_from_string('A_p+_B_[x]')
    c_pplus_b = rfs.reaction_from_string('C_p+_B_[x]')


    rxncon = rxs.RxnConSystem([a_pplus_b, c_pplus_b], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    universal_specA = spe.Specification("A", None, None, None)
    universal_specB = spe.Specification("B", None, None, None)
    universal_specC = spe.Specification("C", None, None, None)

    specificationBp = spe.Specification("B", None, None, 'x')

    expected_mol_def_B = mol.MoleculeDefinition(universal_specB, {mol.ModificationDefinition(specificationBp, {mol.Modifier.unmodified, mol.Modifier.phosphorylated})},
                                                set(), mol.LocalizationDefinition(set()))

    assert not mol_defs.molecule_definition_for_name("A").association_defs
    assert not mol_defs.molecule_definition_for_name("A").modification_defs

    assert not mol_defs.molecule_definition_for_name("C").association_defs
    assert not mol_defs.molecule_definition_for_name("C").modification_defs

    assert not mol_defs.molecule_definition_for_name("B").association_defs
    assert mol_defs.molecule_definition_for_name("B").localization_def == expected_mol_def_B.localization_def

def test_molecule_defs_from_rxncon_different_modifiation_at_same_residue():
    a_pplus_b = rfs.reaction_from_string('A_p+_B_[(x)]')
    c_pplus_b = rfs.reaction_from_string('C_ub+_B_[(x)]')


    rxncon = rxs.RxnConSystem([a_pplus_b, c_pplus_b], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    universal_specB = spe.Specification("B", None, None, None)

    specificationBp = spe.Specification("B", None, None, 'x')

    expected_mol_def_B = mol.MoleculeDefinition(universal_specB, {mol.ModificationDefinition(specificationBp,
                                                                                             {mol.Modifier.unmodified,
                                                                                              mol.Modifier.phosphorylated,
                                                                                              mol.Modifier.ubiquitinated})},
                                                set(), mol.LocalizationDefinition(set()))

    assert not mol_defs.molecule_definition_for_name("A").association_defs
    assert not mol_defs.molecule_definition_for_name("A").modification_defs

    assert not mol_defs.molecule_definition_for_name("C").association_defs
    assert not mol_defs.molecule_definition_for_name("C").modification_defs

    assert not mol_defs.molecule_definition_for_name("B").association_defs
    assert mol_defs.molecule_definition_for_name("B").localization_def == expected_mol_def_B.localization_def
    assert list(mol_defs.molecule_definition_for_name("B").modification_defs)[0].spec == list(expected_mol_def_B.modification_defs)[0].spec
    assert len(list(mol_defs.molecule_definition_for_name("B").modification_defs)[0].valid_modifiers.intersection(list(expected_mol_def_B.modification_defs)[0].valid_modifiers)) == 3

def test_molecule_defs_from_rxncon_binding_same_domain():
    a_pplus_b = rfs.reaction_from_string('A_ppi_B_[x]')
    c_pplus_b = rfs.reaction_from_string('C_ppi_B_[x]')


    rxncon = rxs.RxnConSystem([a_pplus_b, c_pplus_b], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    universal_specB = spe.Specification("B", None, None, None)

    specificationBbound = spe.Specification("B", "x", None, None)
    specificationAiB = spe.Specification("A", "Bassoc", None, None)
    specificationCiB = spe.Specification("C", "Bassoc", None, None)

    expected_mol_def_B = mol.MoleculeDefinition(universal_specB, set(),
                                                {mol.AssociationDefinition(specificationBbound, {specificationAiB, specificationCiB})}, mol.LocalizationDefinition(set()))


    assert len(mol_defs.molecule_definition_for_name("B").association_defs) == 1
    assert list(mol_defs.molecule_definition_for_name("B").association_defs)[0].spec == list(expected_mol_def_B.association_defs)[0].spec
    assert len(list(mol_defs.molecule_definition_for_name("B").association_defs)[0].valid_partners) == 2
    assert len(list(mol_defs.molecule_definition_for_name("B").association_defs)[0].valid_partners.intersection(list(expected_mol_def_B.association_defs)[0].valid_partners)) == 2

