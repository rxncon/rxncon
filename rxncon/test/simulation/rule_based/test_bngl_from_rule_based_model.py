from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon
from rxncon.simulation.rule_based.bngl_from_rule_based_model import bngl_from_rule_based_model


def test_simple_system() -> None:
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                             A_[b]_ppi-_B_[a]
                                             C_p+_A_[(r)]; k+ A_[b]--B_[a]
                                             D_p-_A_[(r)]""").rxncon_system)

    expected_bngl = '''begin model
begin parameters
NumA      1000
NumB      1000
NumC      1000
NumD      1000
k_1       1.0		#  A_[b]_ppi+_B_[a]
k_2       1.0		#  A_[b]_ppi-_B_[a]
k_3_1     1.0		#  C_p+_A_[(r)]
k_3_2     1.0		#  C_p+_A_[(r)]
k_4       1.0		#  D_p-_A_[(r)]
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
# Rule 1. rxn: A_[b]_ppi+_B_[a], quant_cont: UniversalSet
A(bD,rR~p) + B(aD) -> A(bD!1,rR~p).B(aD!1)   k_1

# Rule 2. rxn: A_[b]_ppi-_B_[a], quant_cont: UniversalSet
A(bD!1).B(aD!1) -> A(bD) + B(aD)   k_2

# Rule 3. rxn: C_p+_A_[(r)], quant_cont: A@1_[b]--B@2_[a]
A(bD!1,rR~0).B(aD!1) + C() -> A(bD!1,rR~p).B(aD!1) + C()   k_3_1

# Rule 4. rxn: C_p+_A_[(r)], quant_cont: !(A@1_[b]--B@2_[a])
A(bD,rR~0) + C() -> A(bD,rR~p) + C()   k_3_2

# Rule 5. rxn: D_p-_A_[(r)], quant_cont: UniversalSet
A(rR~p) + D() -> A(rR~0) + D()   k_4

end reaction rules

end model

simulate_nf({t_end=>100,n_steps=>10});
'''
    assert bngl_from_rule_based_model(rbm) == expected_bngl
