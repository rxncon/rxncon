import rxncon.input.quick.quick as qui

def test_quick_single_reaction():
    quick = qui.Quick("A_ppi_B")

    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 1
    assert len(quick._contingencies) == 0

def test_quick_single_reaction_simple_cont():
    quick = qui.Quick("A_ppi_B; ! A--C")
    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 1
    assert len(quick._contingencies) == 1

def test_quick_single_reaction_simple_bool():
    quick = qui.Quick("""A_ppi_B; ! <bool>
                        <bool>; AND A--C
                        <bool>; AND A--D
                        <bool>; AND B--C""")
    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 1
    assert len(quick._contingencies) == 1

def test_quick_multiple_reactions_contingencies():
    quick = qui.Quick("""A_ppi_B; ! A--C
                         B_p+_C; ! <bool>
                         <bool>; OR B--D
                         <bool>; OR <bool1>
                         <bool1>; AND B--E
                         <bool1>; AND B--F
                        """)

    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 2
    assert len(quick._contingencies) == 2

def test():
    quick = qui.Quick("A_ppi_B; ! A_[n]--[m]")
    pass