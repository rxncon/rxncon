import rxncon.core.effector as eff
import rxncon.syntax.rxncon_from_string as rfs


def test_effector_states_property():
    state_a1 = rfs.state_from_string('A--C')
    state_a2 = rfs.state_from_string('A-{p}')
    state_b1 = rfs.state_from_string('B-{ub}')
    state_b2 = rfs.state_from_string('B--D')

    effector = eff.OrEffector(eff.AndEffector(eff.StateEffector(state_a1),
                                              eff.StateEffector(state_a2)),
                              eff.AndEffector(eff.StateEffector(state_b1),
                                              eff.StateEffector(state_b2)))

    assert all(x in effector.states for x in [state_a1, state_a2, state_b1, state_b2])
    assert all(x in [state_a1, state_a2, state_b1, state_b2] for x in effector.states)
