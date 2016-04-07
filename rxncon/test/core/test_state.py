import rxncon.core.state as sta
import rxncon.syntax.rxncon_from_string as rfs


def test_superspec_subspec_covalent_modification():
    assert rfs.state_from_string('A-{p}').is_superspecification_of(rfs.state_from_string('A_[n]-{p}'))
    assert not rfs.state_from_string('A-{p}').is_subspecification_of(rfs.state_from_string('A_[n]-{p}'))

    assert rfs.state_from_string('A_[n]-{p}').is_subspecification_of(rfs.state_from_string('A-{p}'))
    assert not rfs.state_from_string('A_[n]-{p}').is_superspecification_of(rfs.state_from_string('A-{p}'))


def test_superspec_subspec_ppi():
    # Superspecification
    assert rfs.state_from_string('A--B').is_superspecification_of(rfs.state_from_string('A_[n]--B'))
    assert rfs.state_from_string('A--B').is_superspecification_of(rfs.state_from_string('A--B_[m]'))
    assert rfs.state_from_string('A--B').is_superspecification_of(rfs.state_from_string('A_[n]--B_[m]'))

    assert not rfs.state_from_string('A_[n]--B').is_superspecification_of(rfs.state_from_string('A--B_[m]'))
    assert not rfs.state_from_string('A--B_[m]').is_superspecification_of(rfs.state_from_string('A_[n]--B'))
    assert rfs.state_from_string('A_[n]--B').is_superspecification_of(rfs.state_from_string('A_[n]--B_[m]'))
    assert rfs.state_from_string('A--B_[m]').is_superspecification_of(rfs.state_from_string('A_[n]--B_[m]'))

    # Subspecification
    assert rfs.state_from_string('A_[n]--B').is_subspecification_of(rfs.state_from_string('A--B'))
    assert rfs.state_from_string('A--B_[m]').is_subspecification_of(rfs.state_from_string('A--B'))
    assert rfs.state_from_string('A_[n]--B_[m]').is_subspecification_of(rfs.state_from_string('A--B'))

    assert not rfs.state_from_string('A--B_[m]').is_subspecification_of(rfs.state_from_string('A_[n]--B'))
    assert not rfs.state_from_string('A_[n]--B').is_subspecification_of(rfs.state_from_string('A--B_[m]'))
    assert rfs.state_from_string('A_[n]--B_[m]').is_subspecification_of(rfs.state_from_string('A_[n]--B'))
    assert rfs.state_from_string('A_[n]--B_[m]').is_subspecification_of(rfs.state_from_string('A--B_[m]'))

def test_superspec_subspec_component_state():
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A'))

    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A--B'))
    assert rfs.state_from_string('B').is_superspecification_of(rfs.state_from_string('A--B'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d]--B'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s]--B'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s(r)]--B'))

    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A-{p}'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d]-{p}'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s]-{p}'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s(r)]-{p}'))
    assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[(r)]-{p}'))

    assert rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A--B'))
    assert not rfs.state_from_string('B').is_subspecification_of(rfs.state_from_string('A--B'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d]--B'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s]--B'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s(r)]--B'))

    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A-{p}'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d]-{p}'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s]-{p}'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s(r)]-{p}'))
    assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[(r)]-{p}'))
