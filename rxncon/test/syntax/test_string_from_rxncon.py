import rxncon.syntax.rxncon_from_string as rfs


reaction_strings = [
    'Kss1_[CD]_ppi_Ste7_[BD:MAPK]',
    'Ste2_[CyT]_ppi_Sst2_[n/DEP]',
    'Ste7_[KD]_P+_Kss1_[(T183)]',
    'Pkh1_[KD]_P+_Pkc1_[AL(T983)]',
    'Msg5_[PD]_P-_Fus3_[(Y182)]',
    'Sac7_[GAP]_GAP_Rho1_[GnP]',
    'Hog1_[KD]_P+_Sko1_[n(S126)]'
]


state_strings = [
    'MFalpha--Ste2_[Receptor]',
    'Sst2_[(S539)]-{P}',
    'Pkc1_[AL(T983)]-{P}',
    'Rom2_[n]--Slg1_[CyT]',
    'Hog1_[PBD-2]--Pbs2_[HBD-1]',
    'Opy2_[BD:Ste50]--Ste50_[RA]',
    'Dig1_[Dsite]--Kss1_[CD/7m]',
    'Hot1_[(S360)]-{P}'
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


