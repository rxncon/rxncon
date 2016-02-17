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


# TESTING MOLECULE DEFINITION CREATION BY MOLECULEDEFINITIONSUPERVISOR
def test_molecule_definitions_no_contingencies():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c, b_ppi_e, b_pplus_e], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

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

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

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

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

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

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

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

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

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


# TESTING EFFECTOR TO STATES
def test_set_of_states_from_effector_FOR_state_effector():
    a_ppi_c = rfs.reaction_from_string("A_ppi_C")
    a_dash_d = rfs.state_from_string('A--D')

    cont = con.Contingency(a_ppi_c, con.ContingencyType.requirement, eff.StateEffector(a_dash_d))

    set_of_state_effector_a_dash_d = mfr.set_of_states_from_effector(cont.effector)
    assert set_of_state_effector_a_dash_d.is_equivalent_to(venn.PropertySet(a_dash_d))


def test_set_of_states_from_effector_FOR_and_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector)

    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_effector_FOR_or_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR A--C
                        <comp>; OR C--E
                        <comp>; OR B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_effector_FOR_and_or_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND <c1>
                        <comp>; AND <c2>
                        <c1>; OR A--C
                        <c1>; OR C--E
                        <c2>; OR B--F
                        <c2>; OR B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Union(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_effector_FOR_or_and_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; AND B--F
                        <c2>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Intersection(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_effector_FOR_or_and_not_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; NOT B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_d = rfs.state_from_string("B--D")


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Complement(venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_effector_FOR_or_and_not_complex_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; NOT <c3>
                        <c3>; AND B--F
                        <c3>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = mfr.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Complement(venn.Intersection(venn.PropertySet(expected_b_dash_d),
                                                                                                          venn.PropertySet(expected_b_dash_f)))))

# TESTING CONTINGENCIES TO SETS OF STATES
def test_set_of_states_FOR_contingencies_strict():
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

    strict_contingencies_state_set_b_ppi_e = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    strict_contingencies_state_set_a_ppi_b = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[1]))

    strict_contingencies_state_set_a_ppi_c = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[2]))

    strict_contingencies_state_set_b_pplus_e = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[3]))

    expected_b_ppi_e_strict_cont = venn.Intersection(venn.PropertySet(e_pplus), venn.PropertySet(a_dash_b))
    assert strict_contingencies_state_set_b_ppi_e.is_equivalent_to(expected_b_ppi_e_strict_cont)

    expected_a_ppi_b_strict_cont = venn.Intersection(venn.UniversalSet(),venn.Complement(venn.PropertySet(a_dash_c)))
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(expected_a_ppi_b_strict_cont)
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(venn.Complement(venn.PropertySet(a_dash_c)))

    expected_a_ppi_c_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_a_ppi_c.is_equivalent_to(expected_a_ppi_c_strict_cont)

    expected_b_pplus_e_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_b_pplus_e.is_equivalent_to(expected_b_pplus_e_strict_cont)


def test_set_of_states_from_contingencies_FOR_quant():
    # todo
    pass


def test_state_set_from_contingencies_FOR_AND_complex():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")


    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    assert strict_cont_state_set.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c), venn.PropertySet(expected_c_dash_e)),
                                                                    venn.PropertySet(expected_b_dash_f)))

def test_source_set_of_states_FOR_reaction():
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

    # todo: B_pt_E are two reactions in one B_p+_E -> E_[Bside] and E_p-_B -> B_[Eside]
    # todo: B_[n]_apt_B_[m] auto phosphortransfer B is the same molecule B_[n]_p+_B_[m] -> B_[m] and B_[m]_p-_B_[n] -> B_B[n]
    # todo: B_apt_B auto phosphortransfer B is the same molecule B_p+_B -> B_[Bsite1] and B_p-_B -> B_B[Site2]

    assert set_a_ppi_b.is_equivalent_to(venn.Complement(venn.PropertySet(a_dash_b)))
    assert set_b_pplus_e.is_equivalent_to(venn.Complement(venn.PropertySet(e_pplus)))
    assert set_b_pminus_e.is_equivalent_to(venn.PropertySet(e_pplus))

    #very nice that B should be phosphorilated even if there is no reaction for it (later we should test this)
    assert set_b_pt_e.is_equivalent_to(venn.Intersection(venn.Complement(venn.PropertySet(e_pplus)), venn.PropertySet(b_pplus)))


def test_set_of_instances_from_molecule_def_and_set_of_states_FOR_ppi_no_contingency():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    rxncon = rxs.RxnConSystem([a_ppi_b], [])

    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_instances_set = mfr.set_of_instances_from_molecule_def_and_set_of_states(mol_defs.molecule_definition_for_name("A"), strict_cont_state_set)

    assert strict_instances_set.is_equivalent_to(venn.UniversalSet())


def test_set_of_instances_FROM_molecule_def_and_set_of_states_FOR_requirement_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_instances_set = mfr.set_of_instances_from_molecule_def_and_set_of_states(mol_defs.molecule_definition_for_name("A"), strict_cont_state_set)

    assoc_defs = [assoc_def for assoc_def in mol_defs.molecule_definition_for_name("A").association_defs if assoc_def.spec.domain == "Cassoc"]
    assoc_def = assoc_defs[0]

    assert isinstance(assoc_def, mol.AssociationDefinition)
    assert assoc_def.spec == spe.Specification('A', 'Cassoc', None, None)

    expected_assoc_instance = mol.AssociationInstance(assoc_def,
                                                      mol.OccupationStatus.occupied_known_partner,
                                                      spe.Specification('C', 'Aassoc', None, None))

    assert strict_instances_set.is_equivalent_to(venn.PropertySet(expected_assoc_instance))


def test_set_of_instances_from_molecule_def_and_set_of_states_FOR_inhibition_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = mfr.MoleculeDefinitionSupervisor(rxncon)

    strict_cont_state_set = mfr.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_instances_set = mfr.set_of_instances_from_molecule_def_and_set_of_states(mol_defs.molecule_definition_for_name("A"), strict_cont_state_set)

    assoc_defs = [assoc_def for assoc_def in mol_defs.molecule_definition_for_name("A").association_defs if assoc_def.spec.domain == "Cassoc"]
    assoc_def = assoc_defs[0]
    assert isinstance(assoc_def, mol.AssociationDefinition)
    assert assoc_def.spec == spe.Specification('A', 'Cassoc', None, None)

    expected_assoc_instance = mol.AssociationInstance(assoc_def,
                                                      mol.OccupationStatus.not_occupied,
                                                      spe.Specification('C', 'Aassoc', None, None))

    assert strict_instances_set.is_equivalent_to(venn.PropertySet(expected_assoc_instance))
