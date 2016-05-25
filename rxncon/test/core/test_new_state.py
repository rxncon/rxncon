import rxncon.core.state as sta

def test_state():
    sta.state_from_string('[INPUT]')
    state = sta.state_from_string('A_[m]--[n]')
    str(state)
    sta.state_from_string('[INPUT]')
    sta.state_from_string('A--B_[m]')
    state1 = sta.state_from_string('A-{P}')
    state2 = sta.state_from_string('A')
    state3 = sta.state_from_string('A--B')
    state4 = sta.state_from_string('A_[d]--B')

    state2.is_superspecification_of(state1)
    state3.is_superspecification_of(state4)
    state4.is_superspecification_of(state3)
    state2.is_superspecification_of(state3)
    state1.is_superspecification_of(state2)
    state1.is_subspecification_of(state1)

    sta.state_from_string('A--B').is_superspecification_of(sta.state_from_string('A_[n]--B'))
    sta.state_from_string('A--B').is_superspecification_of(sta.state_from_string('A--B_[m]'))
    sta.state_from_string('A--B').is_superspecification_of(sta.state_from_string('A_[n]--B_[m]'))

    sta.state_from_string('A_[n]--B_[m]').is_subspecification_of(sta.state_from_string('A--B'))
    sta.state_from_string('A_[n]--B').is_subspecification_of(sta.state_from_string('A--B'))
    sta.state_from_string('A--B_[m]').is_subspecification_of(sta.state_from_string('A--B'))
    pass