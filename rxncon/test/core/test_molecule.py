from rxncon.input.quick.quick import Quick
from rxncon.semantics.molecule import mol_def_for_component, Mol, InconsistentStructureError, MutualExclusivityError
from rxncon.core.spec import mol_spec_from_string
from rxncon.core.state import state_from_string
from rxncon.util.utils import elems_eq

import pytest

def test_mol_defs():
    rxncon_sys = Quick('''A_trsl_BmRNA
                          C_p+_B@5_[(r1)]
                          D_p+_B_[(r2)] ; ! B@4_[(r1)]-{p}
                          D_[x]_ppi_B_[y]
                          E_ub+_B_[(r1)]''').rxncon_system

    moldef = mol_def_for_component(rxncon_sys, mol_spec_from_string('B'))

    assert elems_eq(list(moldef.valid_states_by_spec(mol_spec_from_string('B_[(r1)]')).values())[0],
                    [state_from_string('B_[(r1)]-{0}'), state_from_string('B_[(r1)]-{p}'),
                     state_from_string('B_[(r1)]-{ub}')])

    mol = Mol(moldef, {})
    assert mol.component == mol_spec_from_string('B')
    mol.add_state(mol_spec_from_string('B_[(r1)]'), state_from_string('B@4_[(r1)]-{p}'))
    assert mol.component == mol_spec_from_string('B@4')

    with pytest.raises(InconsistentStructureError):
        mol.add_state(mol_spec_from_string('B_[y]'), state_from_string('D_[x]--B@7_[y]'))

    with pytest.raises(MutualExclusivityError):
        mol.add_state(mol_spec_from_string('B_[(r1)]'), state_from_string('B_[(r1)]-{ub}'))
