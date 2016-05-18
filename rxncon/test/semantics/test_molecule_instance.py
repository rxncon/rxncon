import pytest
from collections import namedtuple

from rxncon.simulation.rule_based.molecule_from_string import mol_instance_from_string

MoleculeInstanceSortingTestCase = namedtuple('MoleculeInstanceSortingTestCase', ["actual_molecule_instances",
                                                                                 "expected_molecule_instances_sorting"])


@pytest.fixture
def the_case_molecule_instances():
    return [MoleculeInstanceSortingTestCase([mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:'),
                                             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
                                             mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:')],

                                            [mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:'),
                                             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
                                             mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:')]
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:'),
                                             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
                                             mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:')],

                                            [mol_instance_from_string('A#ass/A_[dB]:B_[dA]', 'A#ass/A_[dB]:'),
                                             mol_instance_from_string('B#ass/B_[dA]:A_[dB]', 'B#ass/B_[dA]:'),
                                             mol_instance_from_string('C#ass/C_[dA]:A_[dC]', 'C#ass/C_[dA]:')],
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#ass/A_[dE]:E_[dA]', 'A#ass/A_[dE]:'),
                                             mol_instance_from_string('A#ass/A_[dD]:D_[dA]', 'A#ass/A_[dD]:'),
                                             mol_instance_from_string('A#ass/A_[dC]:C_[dA]', 'A#ass/A_[dC]:')],

                                            [mol_instance_from_string('A#ass/A_[dC]:C_[dA]', 'A#ass/A_[dC]:'),
                                             mol_instance_from_string('A#ass/A_[dD]:D_[dA]', 'A#ass/A_[dD]:'),
                                             mol_instance_from_string('A#ass/A_[dE]:E_[dA]', 'A#ass/A_[dE]:')],
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#ass/A_[dE/sE]:E_[dA]', 'A#ass/A_[dE/sE]:'),
                                             mol_instance_from_string('A#ass/A_[dD/sD]:D_[dA]', 'A#ass/A_[dD/sD]:'),
                                             mol_instance_from_string('A#ass/A_[dC/sC]:C_[dA]', 'A#ass/A_[dC/sC]:')],

                                            [mol_instance_from_string('A#ass/A_[dC/sC]:C_[dA]', 'A#ass/A_[dC/sC]:'),
                                             mol_instance_from_string('A#ass/A_[dD/sD]:D_[dA]', 'A#ass/A_[dD/sD]:'),
                                             mol_instance_from_string('A#ass/A_[dE/sE]:E_[dA]', 'A#ass/A_[dE/sE]:')],
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:'),
                                             mol_instance_from_string('A#ass/A_[rD]:D_[dA]', 'A#ass/A_[rD]:'),
                                             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:')],

                                            [mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:'),
                                             mol_instance_from_string('A#ass/A_[rD]:D_[dA]', 'A#ass/A_[rD]:'),
                                             mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:')],
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:'),
                                             mol_instance_from_string('A#', 'A#'),
                                             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:')],

                                            [mol_instance_from_string('A#', 'A#'),
                                             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:'),
                                             mol_instance_from_string('A#ass/A_[rE]:E_[dA]', 'A#ass/A_[rE]:')],
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#ass/A_[rE]:E_[dA],mod/A_[(r)]:u~p', 'A#ass/A_[rE]:,mod/A_[(r)]:u'),
                                             mol_instance_from_string('A#', 'A#'),
                                             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:')],

                                            [mol_instance_from_string('A#', 'A#'),
                                             mol_instance_from_string('A#ass/A_[rC]:C_[dA]', 'A#ass/A_[rC]:'),
                                             mol_instance_from_string('A#ass/A_[rE]:E_[dA],mod/A_[(r)]:u~p',
                                                                      'A#ass/A_[rE]:,mod/A_[(r)]:u')],
                                            ),

            MoleculeInstanceSortingTestCase([mol_instance_from_string('A#mod/A_[(r)]:u~p', 'A#mod/A_[(r)]:p'),
                                             mol_instance_from_string('A#mod/A_[(r)]:u~ub', 'A#mod/A_[(r)]:ub'),
                                             mol_instance_from_string('A#mod/A_[(r)]:u~ub~p', 'A#mod/A_[(r)]:u')],

                                            [mol_instance_from_string('A#mod/A_[(r)]:u~p', 'A#mod/A_[(r)]:p'),
                                             mol_instance_from_string('A#mod/A_[(r)]:u~ub~p', 'A#mod/A_[(r)]:u'),
                                             mol_instance_from_string('A#mod/A_[(r)]:u~ub', 'A#mod/A_[(r)]:ub')]
                                            )

            ]


def test_molecule_instance_sorting(the_case_molecule_instances,):
    for the_case in the_case_molecule_instances:
        assert sorted(the_case.actual_molecule_instances) == the_case.expected_molecule_instances_sorting
