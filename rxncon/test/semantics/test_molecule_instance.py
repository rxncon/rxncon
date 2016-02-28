import pytest
import rxncon.venntastic.sets as venn
import rxncon.core.rxncon_system as rxs
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.semantics.molecule_instance as mins
import rxncon.semantics.molecule_definition as mdef
import rxncon.semantics.molecule_instance_from_rxncon as mifr
import rxncon.semantics.molecule_definition_from_rxncon as mdfr


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
