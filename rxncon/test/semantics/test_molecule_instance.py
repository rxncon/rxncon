import pytest
import rxncon.venntastic.sets as venn
import rxncon.core.rxncon_system as rxs
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.semantics.molecule_instance as mins
import rxncon.semantics.molecule_definition as mdef
import rxncon.semantics.molecule_instance_from_rxncon as mifr
import rxncon.semantics.molecule_definition_from_rxncon as mdfr
import rxncon.semantics.molecule_definition as mdf
import rxncon.core.specification as spe


def test_commutative_diagram_for_complement_operation(mol_def, state_sets):
    for state_set in state_sets:
        mol_props_then_complement = venn.Complement(mifr.property_set_from_mol_def_and_state_set(mol_def, state_set))
        complement_then_mol_props = mifr.property_set_from_mol_def_and_state_set(mol_def, venn.Complement(state_set))

        print()
        print(mol_props_then_complement.simplified_form())
        print()
        print(complement_then_mol_props.simplified_form())

        #assert mol_props_then_complement.is_equivalent_to(complement_then_mol_props)


@pytest.fixture
def mol_def() -> mdef.MoleculeDefinition:
    reactions = [
        rfs.reaction_from_string('B_ppi_A_[x]'),
        rfs.reaction_from_string('C_ppi_A_[x]')
    ]

    rxnsys = rxs.RxnConSystem(reactions, [])
    mol_def_supervisor = mdfr.MoleculeDefinitionSupervisor(rxnsys)

    return mol_def_supervisor.mol_def_for_name('A')


@pytest.fixture
def state_sets():
    return [
        venn.PropertySet(rfs.state_from_string('A--B'))
    ]

@pytest.fixture
def specification_generic_domain_B():
    return spe.Specification("B", 'generic', None, None)

@pytest.fixture
def spec_A_dB():
    return spe.Specification("A", 'dB', None, None)

@pytest.fixture
def spec_A_dC():
    return spe.Specification("A", 'dC', None, None)

@pytest.fixture
def spec_A_sB():
    return spe.Specification("A", None, 'sB', None)

@pytest.fixture
def spec_A_sC():
    return spe.Specification("A", None, 'sC', None)

@pytest.fixture
def spec_A_rB():
    return spe.Specification("A", None, None, 'rB')

@pytest.fixture
def spec_A_rC():
    return spe.Specification("A", None, None, 'rC')

@pytest.fixture
def spec_A_dB_sC():
    return spe.Specification("A", "dB", "sC", None)

@pytest.fixture
def spec_A_dB_sD():
    return spe.Specification("A", "dB", "sD", None)


@pytest.fixture
def spec_A_dB_rC():
    return spe.Specification("A", "dB", None, "rC")


@pytest.fixture
def spec_A_dB_rD():
    return spe.Specification("A", "dB", None, "sD")


@pytest.fixture
def spec_A_sB_rC():
    return spe.Specification("A", None, "sB", "rC")


@pytest.fixture
def spec_A_sB_rD():
    return spe.Specification("A", None, "sB", "rD")


def assoc_instance_A(spec_A_1, spec_A_2, specification_generic_domain_B):
    return {mins.AssociationPropertyInstance(mdf.AssociationPropertyDefinition(spec_A_1, {specification_generic_domain_B}),
                                             mdf.OccupationStatus.occupied_known_partner,
                                             specification_generic_domain_B),
            mins.AssociationPropertyInstance(mdf.AssociationPropertyDefinition(spec_A_2, {specification_generic_domain_B}),
                                             mdf.OccupationStatus.occupied_known_partner,
                                             specification_generic_domain_B)}


def expected_ordering_assoc_instance(spec_A_1, spec_A_2, specification_generic_domain_B):
    return [mins.AssociationPropertyInstance(mdf.AssociationPropertyDefinition(spec_A_1, {specification_generic_domain_B}),
                                             mdf.OccupationStatus.occupied_known_partner,
                                             specification_generic_domain_B),
            mins.AssociationPropertyInstance(mdf.AssociationPropertyDefinition(spec_A_2, {specification_generic_domain_B}),
                                             mdf.OccupationStatus.occupied_known_partner,
                                             specification_generic_domain_B)
            ]

@pytest.fixture
def association_property_instances(spec_A_dB, spec_A_dC, spec_A_sB, spec_A_sC, spec_A_rB, spec_A_rC,
                                   spec_A_dB_sC, spec_A_dB_sD, spec_A_dB_rC, spec_A_dB_rD, spec_A_sB_rC, spec_A_sB_rD,
                                   specification_generic_domain_B):
    return [assoc_instance_A(spec_A_dC, spec_A_dB, specification_generic_domain_B),
            assoc_instance_A(spec_A_sC, spec_A_sB, specification_generic_domain_B),
            assoc_instance_A(spec_A_rC, spec_A_rB, specification_generic_domain_B),
            assoc_instance_A(spec_A_dB_sD, spec_A_dB_sC, specification_generic_domain_B),
            assoc_instance_A(spec_A_dB_rD, spec_A_dB_rC, specification_generic_domain_B),
            assoc_instance_A(spec_A_sB_rD, spec_A_sB_rC, specification_generic_domain_B)]

@pytest.fixture
def expected_association_domain_ordering(spec_A_dB, spec_A_dC, spec_A_sB, spec_A_sC, spec_A_rB, spec_A_rC,
                                         spec_A_dB_sC, spec_A_dB_sD, spec_A_dB_rC, spec_A_dB_rD, spec_A_sB_rC, spec_A_sB_rD,
                                         specification_generic_domain_B):
    return [expected_ordering_assoc_instance(spec_A_dB, spec_A_dC, specification_generic_domain_B),
            expected_ordering_assoc_instance(spec_A_sB, spec_A_sC, specification_generic_domain_B),
            expected_ordering_assoc_instance(spec_A_rB, spec_A_rC, specification_generic_domain_B),
            expected_ordering_assoc_instance(spec_A_dB_sC, spec_A_dB_sD, specification_generic_domain_B),
            expected_ordering_assoc_instance(spec_A_dB_rC, spec_A_dB_rD, specification_generic_domain_B),
            expected_ordering_assoc_instance(spec_A_sB_rC, spec_A_sB_rD, specification_generic_domain_B)]

def test_association_property_instance_ordering(association_property_instances, expected_association_domain_ordering):
    assert all(sorted(instance) == expected_association_domain_ordering[i] for i, instance in enumerate(association_property_instances) )


# def molecule_definition_A():
#     spec_A = spe.Specification("A", None, None, None)
#     mol_def_A = mdf.MoleculeDefinition(
#         spec_A,
#         {mdf.ModificationPropertyDefinition(spec_A_psite, {mdf.Modifier.unmodified, mdf.Modifier.phosphorylated})},
#         {mdf.AssociationPropertyDefinition(spec_A_Bassoc, {spec_B_Aassoc})},
#         None
#     )
def modification_definition_A():
    spec_A_psite = spe.Specification('A', None, None, 'psite')

    spec_B_Aassoc = spe.Specification("B", 'Aassoc', None, None)

    spec_A_psite = spe.Specification('A', None, None, 'psite')



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
