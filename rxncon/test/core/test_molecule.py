from rxncon.input.quick.quick import Quick
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.semantics.molecule import MolDef, mol_def_for_component
from rxncon.core.spec import mol_spec_from_string


def test_mol_defs():
    rxncon_sys = Quick('''A_trsl_BmRNA
                          C_p+_B_[(r1)]
                          D_p+_B_[(r2)] ; ! B_[(r1)]-{p}
                          D_[x]_ppi_B_[y]
                          E_ub+_B_[(r1)]''').rxncon_system

    print(mol_def_for_component(rxncon_sys, mol_spec_from_string('B')))
