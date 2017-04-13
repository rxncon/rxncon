import pytest
from typing import Dict
from rxncon.core.reaction import reaction_from_str, initialize_reaction_defs
from rxncon.core.state import state_from_str
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.simulation.boolean.boolean_model import boolean_model_from_rxncon, SmoothingStrategy
from rxncon.test.simulation.boolean.utils import target_from_str
from rxncon.venntastic.sets import UniversalSet, EmptySet


def check_motif(motif: Dict[str, bool], in_state: Dict[str, bool], out_state: Dict[str, bool],
                smoothing: SmoothingStrategy=SmoothingStrategy.smooth_production_sources) -> bool:
    rxncon_sys = RxnConSystem([reaction_from_str(x) for x, _ in motif.items()], [])
    boolean_model = boolean_model_from_rxncon(rxncon_sys, smoothing)

    for rxn, val in motif.items():
        if val:
            boolean_model.update_rule_by_target(target_from_str(rxn)).factor = UniversalSet()
        else:
            boolean_model.update_rule_by_target(target_from_str(rxn)).factor = EmptySet()

    for state, val in in_state.items():
        boolean_model.initial_conditions.set_target(target_from_str(state), val)

    steady_state = boolean_model.calc_steady_state()

    for state, val in out_state.items():
        if steady_state[target_from_str(state)] != val:
            print()
            print('in  {}'.format(in_state))
            print('out {}'.format(out_state))
            print('ss  {}'.format(steady_state.target_to_value))
            print()
            return False

    return True


def test_no_steady_state():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': False
    }

    with pytest.raises(AssertionError):
        assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False},
                           SmoothingStrategy.no_smoothing)

    with pytest.raises(AssertionError):
        assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True},
                           SmoothingStrategy.no_smoothing)


def register_targeted_degradation():
    tdeg_def = {
        '!UID:Reaction': 'targeted-deg',
        '!UID:ReactionKey': 'tdeg',
        '!BidirectionalVerb': 'no',
        '!MolTypeX': 'Protein',
        '!ResolutionX': 'component',
        '!MolTypeY': 'Protein',
        '!ResolutionY': 'residue',
        '!SkeletonRule': '$x%# + $y%#$y%-{p} -> $x%#'
    }

    initialize_reaction_defs([tdeg_def])


def unregister_targeted_degradation():
    initialize_reaction_defs([])


### TEST HELPERS ###


def test_tdeg_registration():
    with pytest.raises(SyntaxError):
        reaction_from_str('A_tdeg_B_[(r)]')

    register_targeted_degradation()

    assert reaction_from_str('A_tdeg_B_[(r)]').degraded_states == [state_from_str('B_[(r)]-{p}')]

    unregister_targeted_degradation()

    with pytest.raises(SyntaxError):
        reaction_from_str('A_tdeg_B_[(r)]')


### MODIFICATION MOTIFS ###


def test_mod_no_reactions():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_pplus():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_pminus():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_pplus_pminus():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_deg():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_deg_pplus():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_deg_pminus():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_deg_pplus_pminus():
    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})


def test_mod_syn():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})


def test_mod_syn_pplus():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})


def test_mod_syn_pminus():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})


def test_mod_syn_pplus_pminus():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})


def test_mod_syn_deg():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})


def test_mod_syn_deg_pplus():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})


def test_mod_syn_deg_pminus():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})


def test_mod_syn_deg_pplus_pminus():
    motif = {
        'X_syn_A': True,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True})


### MODIFICATION MOTIFS TARGETED DEGRADATION ###


def test_mod_tdeg():
    register_targeted_degradation()

    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': False,
        'W_tdeg_A_[(r)]': True,
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})

    unregister_targeted_degradation()


def test_mod_tdeg_pplus():
    register_targeted_degradation()

    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': False,
        'W_tdeg_A_[(r)]': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})

    unregister_targeted_degradation()


def test_mod_tdeg_pminus():
    register_targeted_degradation()

    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': False,
        'Z_p-_A_[(r)]': True,
        'W_tdeg_A_[(r)]': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})

    unregister_targeted_degradation()


def test_mod_tdeg_pplus_pminus():
    register_targeted_degradation()

    motif = {
        'X_syn_A': False,
        'Y_p+_A_[(r)]': True,
        'Z_p-_A_[(r)]': True,
        'W_tdeg_A_[(r)]': True
    }

    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': True, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': True}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})
    assert check_motif(motif, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False}, {'A_[(r)]-{0}': False, 'A_[(r)]-{p}': False})

    unregister_targeted_degradation()


### INTERACTION MOTIFS ###


def test_int_no_reactions():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_ppiplus():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_ppiminus():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_ppiplus_ppiminus():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_deg():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_deg_ppiplus():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': True
    }

    # assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    # assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    # assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    # assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, smoothing=SmoothingStrategy.no_smoothing)
    # assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    # assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    # assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_deg_ppiminus():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_deg_ppiplus_ppiminus():
    motif = {
        'X_syn_A': False,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_ppiplus():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_ppiminus():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_ppiplus_ppiminus():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': False
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_deg():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_deg_ppiplus():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': False,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_deg_ppiminus():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': False,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})


def test_int_syn_deg_ppiplus_ppiminus():
    motif = {
        'X_syn_A': True,
        'A_[b]_ppi+_B_[a]': True,
        'A_[b]_ppi-_B_[a]': True,
        'W_deg_A': True
    }

    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': True}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': True, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': True, 'A_[b]--B_[a]': True})
    assert check_motif(motif, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
    assert check_motif(motif, {'A_[b]--0': False, 'B_[a]--0': False, 'A_[b]--B_[a]': False}, {'A_[b]--0': True, 'B_[a]--0': False, 'A_[b]--B_[a]': False})
