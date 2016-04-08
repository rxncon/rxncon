import pytest
import rxncon.venntastic.sets as venn
import rxncon.core.rxncon_system as rxs
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.semantics.molecule_definition as mdef

from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, mol_instance_from_string
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys


@pytest.fixture
def mol_def() -> mdef.MoleculeDefinition:
    rxnsys = rxs.RxnConSystem([rfs.reaction_from_string('B_ppi_A_[x]'), rfs.reaction_from_string('C_ppi_A_[x]')], [])
    return mol_defs_from_rxncon_sys(rxnsys)[rfs.specification_from_string('A')]


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
             mol_instance_from_string('A#ass/A_[dC]:C_[dA]', 'A#ass/A_[dC]:')],

            [mol_instance_from_string('A#ass/A_[dE/sE]:E_[dA]', 'A#ass/A_[dE/sE]:'),
             mol_instance_from_string('A#ass/A_[dD/sD]:D_[dA]', 'A#ass/A_[dD/sD]:'),
             mol_instance_from_string('A#ass/A_[dC/sC]:C_[dA]', 'A#ass/A_[dC/sC]:')],

            [mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:'),
             mol_instance_from_string('A#ass/A_[rD]:D_[dA]', 'A#ass/A_[rD]:'),
             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:')],

            [mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:'),
             mol_instance_from_string('A#', 'A#'),
             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:')],

            [mol_instance_from_string('A#ass/A_[rE]:E_[dA],mod/A_[(r)]:u~p', 'A#ass/A_[rE]:,mod/A_[(r)]:u'),
             mol_instance_from_string('A#', 'A#'),
             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:')],

            [mol_instance_from_string('A#mod/A_[(r)]:u~p', 'A#mod/A_[(r)]:p'),
             mol_instance_from_string('A#mod/A_[(r)]:u~ub', 'A#mod/A_[(r)]:ub'),
             mol_instance_from_string('A#mod/A_[(r)]:u~ub~p', 'A#mod/A_[(r)]:u')]

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
              mol_instance_from_string('A#ass/A_[dE]:E_[dA]', 'A#ass/A_[dE]:')],

             [mol_instance_from_string('A#ass/A_[dC/sC]:C_[dA]', 'A#ass/A_[dC/sC]:'),
              mol_instance_from_string('A#ass/A_[dD/sD]:D_[dA]', 'A#ass/A_[dD/sD]:'),
              mol_instance_from_string('A#ass/A_[dE/sE]:E_[dA]', 'A#ass/A_[dE/sE]:')],

             [mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:'),
              mol_instance_from_string('A#ass/A_[rD]:D_[dA]', 'A#ass/A_[rD]:'),
              mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:')],

             [mol_instance_from_string('A#', 'A#'),
              mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:'),
              mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:')],

             [mol_instance_from_string('A#', 'A#'),
              mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:'),
              mol_instance_from_string('A#ass/A_[rE]:E_[dA],mod/A_[(r)]:u~p', 'A#ass/A_[rE]:,mod/A_[(r)]:u')],

             [mol_instance_from_string('A#mod/A_[(r)]:u~p', 'A#mod/A_[(r)]:p'),
              mol_instance_from_string('A#mod/A_[(r)]:u~ub~p', 'A#mod/A_[(r)]:u'),
              mol_instance_from_string('A#mod/A_[(r)]:u~ub', 'A#mod/A_[(r)]:ub'),
             ]
            ]

def test_molecule_instance_sorting(molecule_instances, expected_molecule_instances_ordering):
    for i, instances in enumerate(molecule_instances):
        assert sorted(instances) == expected_molecule_instances_ordering[i]
