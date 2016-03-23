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

from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, mol_instance_from_string


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
def molecule_instances():
    return [[mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:'),
             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
             mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:')],

            [mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:'),
             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
             mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:')],

            [mol_instance_from_string('A#ass/A_[dE]:E_[dA]', 'A#ass/A_[dE]:'),
             mol_instance_from_string('A#ass/A_[dD]:D_[dA]', 'A#ass/A_[dD]:'),
             mol_instance_from_string('A#ass/A_[dC]:C_[dA]', 'A#ass/A_[dC]:')]
            ]

@pytest.fixture
def expected_molecule_instances_ordering():
     return [[mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:'),
             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
             mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:')],

             [mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:'),
              mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
              mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:')],

             [mol_instance_from_string('A#ass/A_[dC]:C_[dA]', 'A#ass/A_[dC]:'),
              mol_instance_from_string('A#ass/A_[dD]:D_[dA]', 'A#ass/A_[dD]:'),
              mol_instance_from_string('A#ass/A_[dE]:E_[dA]', 'A#ass/A_[dE]:')]
            ]

def test_molecule_instance_sorting(molecule_instances, expected_molecule_instances_ordering):
    for i, instances in enumerate(molecule_instances):
        assert sorted(instances) == expected_molecule_instances_ordering[i]

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
