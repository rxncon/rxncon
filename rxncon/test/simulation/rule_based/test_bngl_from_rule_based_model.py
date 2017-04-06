from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon
from rxncon.simulation.rule_based.bngl_from_rule_based_model import bngl_from_rule_based_model

def test_simple_system() -> None:
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                             A_[b]_ppi-_B_[a]
                                             C_p+_A_[(r)]
                                             D_p-_A_[(r)]""").rxncon_system)

    expected_bngl = '''begin model
begin parameters
NumA		1000
NumB		1000
NumC		1000
NumD		1000
k		1.0
end parameters

begin molecule types
A(bD,rR~0~p)
B(aD)
C()
D()
end molecule types

begin seed species
A(bD,rR~0)	NumA
B(aD)	NumB
C()	NumC
D()	NumD
end seed species

begin observables

end observables

begin reaction rules
# A_[b]_ppi+_B_[a]
A(bD,rR~p) + B(aD) -> A(bD!1,rR~p).B(aD!1)   k
# A_[b]_ppi-_B_[a]
A(bD!1).B(aD!1) -> A(bD) + B(aD)   k
# C_p+_A_[(r)]
A(rR~0) + C() -> A(rR~p) + C()   k
# D_p-_A_[(r)]
A(rR~p) + D() -> A(rR~0) + D()   k
end reaction rules

end model

simulate_nf({t_end=>100,n_steps=>10});
'''

    assert bngl_from_rule_based_model(rbm) == expected_bngl
