import pytest
from typing import Dict
from rxncon.core.reaction import reaction_from_str
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

