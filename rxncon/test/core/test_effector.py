import rxncon.core.effector as eff
import rxncon.core.state as sta


def test_effector_states_property() -> None:
    state_a1 = sta.state_from_str('A@0--C@2')
    state_a2 = sta.state_from_str('A@0-{p}')
    state_b1 = sta.state_from_str('B@1-{ub}')
    state_b2 = sta.state_from_str('B@1--D@3')

    effector = eff.OrEffector(eff.AndEffector(eff.StateEffector(state_a1),
                                              eff.StateEffector(state_a2)),
                              eff.AndEffector(eff.StateEffector(state_b1),
                                              eff.StateEffector(state_b2)))

    assert all(x in effector.states for x in [state_a1, state_a2, state_b1, state_b2])
    assert all(x in [state_a1, state_a2, state_b1, state_b2] for x in effector.states)
