import rxncon.syntax.rxncon_from_string as rfs


reaction_strings = [
    'Kss1_[CD]_ppi_Ste7_[BD:MAPK]',
    'Ste2_[CyT]_ppi_Sst2_[n/DEP]',
    'Ste7_[KD]_p+_Kss1_[(T183)]',
    'Pkh1_[KD]_p+_Pkc1_[AL(T983)]',
    'Msg5_[PD]_p-_Fus3_[(Y182)]',
    'Sac7_[GAP]_gap_Rho1_[GnP]',
    'Hog1_[KD]_p+_Sko1_[n(S126)]'
]


state_strings = [
    'MFalpha--Ste2_[Receptor]',
    'Sst2_[(S539)]-{p}',
    'Pkc1_[AL(T983)]-{p}',
    'Rom2_[n]--Slg1_[CyT]',
    'Hog1_[PBD-2]--Pbs2_[HBD-1]',
    'Opy2_[BD:Ste50]--Ste50_[RA]',
    'Dig1_[Dsite]--Kss1_[CD/7m]',
    'Hot1_[(S360)]-{p}'
]


def test_string_from_rxncon_reactions_inverse():
    for reaction_string in reaction_strings:
        reaction = rfs.reaction_from_string(reaction_string)

        assert str(reaction) == reaction_string
        assert rfs.reaction_from_string(str(reaction)) == reaction


def test_string_from_rxncon_states_inverse():
    for state_string in state_strings:
        state = rfs.state_from_string(state_string)

        assert str(state) == state_string
        assert rfs.state_from_string(str(state)) == state


