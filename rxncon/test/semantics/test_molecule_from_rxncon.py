import copy

import rxncon.semantics.molecule_from_rxncon as mfr
import rxncon.semantics.molecule as mol
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.rxncon_system as rxs
import rxncon.core.specification as spe
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.venntastic.sets as venn
import rxncon.input.quick.quick as qui


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
    a_pplus_b = rfs.reaction_from_string('A_p+_B_[(x)]')
    c_pplus_b = rfs.reaction_from_string('C_p+_B_[(x)]')


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
    assert list(mol_defs.molecule_definition_for_name("B").modification_defs)[0].spec == list(expected_mol_def_B.modification_defs)[0].spec
    assert len(list(mol_defs.molecule_definition_for_name("B").modification_defs)[0].valid_modifiers.intersection(list(expected_mol_def_B.modification_defs)[0].valid_modifiers)) == 2


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

# TESTING EFFECTOR TO STATES
def test_set_of_states_from_effector_state_effector():
    a_ppi_c = rfs.reaction_from_string("A_ppi_C")
    a_dash_d = rfs.state_from_string('A--D')

    #todo: the StateEffector gets changed due to the application of set_of_states_from_effector
    #todo: set deepcopy
    cont = con.Contingency(a_ppi_c, con.ContingencyType.requirement, eff.StateEffector(a_dash_d))

    expected_a_dash_d = copy.deepcopy(a_dash_d)
    expected_a_dash_d.first_component.domain = "Dassoc"
    expected_a_dash_d.second_component.domain = "Aassoc"

    set_of_state_effector_a_dash_d = mfr.set_of_states_from_effector(cont.effector, cont.target)
    assert set_of_state_effector_a_dash_d.is_equivalent_to(venn.PropertySet(expected_a_dash_d))


def test_set_of_states_from_effector_AND_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    # don't have to deepcopy because two states are different objects anyway
    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_f.first_component.domain = "Fassoc"
    expected_b_dash_f.second_component.domain = "Bassoc"


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector,rxncon.contingencies[0].target)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_effector_OR_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR A--C
                        <comp>; OR C--E
                        <comp>; OR B--F""")

    # don't have to deepcopy because two states are different objects anyway
    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_f.first_component.domain = "Fassoc"
    expected_b_dash_f.second_component.domain = "Bassoc"

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector, rxncon.contingencies[0].target)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_effector_AND_OR_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND <c1>
                        <comp>; AND <c2>
                        <c1>; OR A--C
                        <c1>; OR C--E
                        <c2>; OR B--F
                        <c2>; OR B--D""")

    # don't have to deepcopy because two states are different objects anyway
    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_f.first_component.domain = "Fassoc"
    expected_b_dash_f.second_component.domain = "Bassoc"

    expected_b_dash_d = rfs.state_from_string("B--D")
    expected_b_dash_d.first_component.domain = "Dassoc"
    expected_b_dash_d.second_component.domain = "Bassoc"

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector, rxncon.contingencies[0].target)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Union(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_effector_OR_AND_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; AND B--F
                        <c2>; AND B--D""")

    # don't have to deepcopy because two states are different objects anyway
    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_f.first_component.domain = "Fassoc"
    expected_b_dash_f.second_component.domain = "Bassoc"

    expected_b_dash_d = rfs.state_from_string("B--D")
    expected_b_dash_d.first_component.domain = "Dassoc"
    expected_b_dash_d.second_component.domain = "Bassoc"


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector, rxncon.contingencies[0].target)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Intersection(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_effector_OR_AND_NOT_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; NOT B--D""")

    # don't have to deepcopy because two states are different objects anyway
    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_d = rfs.state_from_string("B--D")
    expected_b_dash_d.first_component.domain = "Dassoc"
    expected_b_dash_d.second_component.domain = "Bassoc"


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector, rxncon.contingencies[0].target)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Complement(venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_effector_OR_AND_NOT_complex_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; NOT <c3>
                        <c3>; AND B--F
                        <c3>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_f.first_component.domain = "Fassoc"
    expected_b_dash_f.second_component.domain = "Bassoc"

    expected_b_dash_d = rfs.state_from_string("B--D")
    expected_b_dash_d.first_component.domain = "Dassoc"
    expected_b_dash_d.second_component.domain = "Bassoc"


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector, rxncon.contingencies[0].target)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Complement(venn.Intersection(venn.PropertySet(expected_b_dash_d),
                                                                                                          venn.PropertySet(expected_b_dash_f)))))

# TESTING CONTINGENCIES TO SETS OF STATES
def test_set_of_states_from_contingencies_strict():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')

    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string("E-{P}")

    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(e_pplus))  # B_ppi_E; ! E-{P}
    cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C

    rxncon = rxs.RxnConSystem([b_ppi_e, a_ppi_b, a_ppi_c, b_pplus_e], [cont_e_pplus, cont_b_dash_e, cont_a_ppi_b])


    expected_a_dash_b = copy.deepcopy(a_dash_b)
    expected_a_dash_b.first_component.domain = "Bassoc"
    expected_a_dash_b.second_component.domain = "Aassoc"

    expected_a_dash_c = copy.deepcopy(a_dash_c)
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"


    expected_e_pplus = copy.deepcopy(e_pplus)
    expected_e_pplus.substrate.residue = "Bsite"

    strict_contingencies_state_set_b_ppi_e = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    strict_contingencies_state_set_a_ppi_b = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[1]))

    strict_contingencies_state_set_a_ppi_c = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[2]))

    strict_contingencies_state_set_b_pplus_e = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[3]))

    expected_b_ppi_e_strict_cont = venn.Intersection(venn.PropertySet(expected_e_pplus), venn.PropertySet(expected_a_dash_b))
    assert strict_contingencies_state_set_b_ppi_e.is_equivalent_to(expected_b_ppi_e_strict_cont)

    expected_a_ppi_b_strict_cont = venn.Intersection(venn.UniversalSet(),venn.Complement(venn.PropertySet(expected_a_dash_c)))
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(expected_a_ppi_b_strict_cont)
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(venn.Complement(venn.PropertySet(expected_a_dash_c)))

    expected_a_ppi_c_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_a_ppi_c.is_equivalent_to(expected_a_ppi_c_strict_cont)

    expected_b_pplus_e_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_b_pplus_e.is_equivalent_to(expected_b_pplus_e_strict_cont)


def test_set_of_states_from_contingencies_quant():
    # todo
    pass


def test_state_set_from_contingencies_from_AND_complex():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")


    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_a_dash_c.first_component.domain = "Cassoc"
    expected_a_dash_c.second_component.domain = "Aassoc"

    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_c_dash_e.first_component.domain = "Eassoc"
    expected_c_dash_e.second_component.domain = "Cassoc"

    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_f.first_component.domain = "Fassoc"
    expected_b_dash_f.second_component.domain = "Bassoc"

    rxncon = quick.rxncon_system

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    assert strict_cont_state_set.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c), venn.PropertySet(expected_c_dash_e)),
                                                                    venn.PropertySet(expected_b_dash_f)))

def test_source_set_of_states_from_reaction():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')
    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string('E-{P}')
    b_pminus_e = rfs.reaction_from_string('B_p-_E')
    b_pt_e = rfs.reaction_from_string('B_pt_E')
    b_pplus = rfs.state_from_string('B-{P}')


    rxncon = rxs.RxnConSystem([a_ppi_b, b_pplus_e, b_pminus_e, b_pt_e], [])

    set_a_ppi_b = mfr.source_set_of_states_from_reaction(rxncon.reactions[0])
    set_b_pplus_e = mfr.source_set_of_states_from_reaction(rxncon.reactions[1])
    set_b_pminus_e = mfr.source_set_of_states_from_reaction(rxncon.reactions[2])
    set_b_pt_e = mfr.source_set_of_states_from_reaction(rxncon.reactions[3])

    expected_a_dash_b = copy.deepcopy(a_dash_b)
    expected_a_dash_b.first_component.domain = "Bassoc"
    expected_a_dash_b.second_component.domain = "Aassoc"

    expected_b_pplus = copy.deepcopy(b_pplus)

    # todo: B_pt_E are two reactions in one B_p+_E -> E_[Bside] and E_p-_B -> B_[Eside]
    # todo: B_[n]_apt_B_[m] auto phosphortransfer B is the same molecule B_[n]_p+_B_[m] -> B_[m] and B_[m]_p-_B_[n] -> B_B[n]
    # todo: B_apt_B auto phosphortransfer B is the same molecule B_p+_B -> B_[Bsite1] and B_p-_B -> B_B[Site2]

    expected_b_pplus.substrate.residue = "Bsite"

    expected_e_pplus = copy.deepcopy(e_pplus)
    expected_e_pplus.substrate.residue = "Bsite"

    assert set_a_ppi_b.is_equivalent_to(venn.Complement(venn.PropertySet(expected_a_dash_b)))
    assert set_b_pplus_e.is_equivalent_to(venn.Complement(venn.PropertySet(expected_e_pplus)))
    assert set_b_pminus_e.is_equivalent_to(venn.PropertySet(expected_e_pplus))

    #very nice that B should be phosphorilated even if there is no reaction for it (later we should test this)
    assert set_b_pt_e.is_equivalent_to(venn.Intersection(venn.Complement(venn.PropertySet(expected_e_pplus)), venn.PropertySet(expected_b_pplus)))


def test_set_of_instances_from_molecule_def_and_set_of_states_ppi_no_contingency():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    rxncon = rxs.RxnConSystem([a_ppi_b], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    strict_spec_set = mfr.set_of_instances_from_molecule_def_and_set_of_states(mol_defs.molecule_definition_for_name("A"), strict_cont_state_set)

    assert strict_spec_set.is_equivalent_to(venn.UniversalSet())


def test_set_of_instances_from_molecule_def_and_set_of_states_requirement_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set = mfr.set_of_instances_from_molecule_def_and_set_of_states(mol_defs.molecule_definition_for_name("A"), strict_cont_state_set)

    assoc_def = [assoc_def for assoc_def in mol_defs.molecule_definition_for_name("A").association_defs if assoc_def.spec.domain == "Cassoc"]

    assert strict_spec_set.is_equivalent_to(venn.PropertySet(
        mol.AssociationInstance(assoc_def[0], mol.OccupationStatus.occupied_known_partner,list(assoc_def[0].valid_partners)[0])))
