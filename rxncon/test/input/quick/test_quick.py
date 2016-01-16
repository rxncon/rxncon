import rxncon.input.quick.quick as qui

# @basti I kind of like the PEP8 guidelines for code formatting, this states you should have 2 line breaks between
#        function definitions and betwee class definitions etc. You can turn PEP8 checking in Pycharm.

def test_quick_single_reaction():
    quick = qui.Quick("A_ppi_B")

    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 1
    assert len(quick._contingencies) == 0
    # @basti You should check that it's also the right reaction that gets created. Also in the other tests.

def test_quick_single_reaction_simple_cont():
    # @basti 'cont' -> 'contingency'. It's the name of a test, should be as descriptive as possible.
    #        also check whether it's the correct contingency that gets created.
    quick = qui.Quick("A_ppi_B; ! A--C")
    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 1
    assert len(quick._contingencies) == 1

def test_quick_single_reaction_simple_bool():
    # @basti check whether the contingency is the one you want. Also 'bool' -> 'boolean_contingency'
    quick = qui.Quick("""A_ppi_B; ! <bool>
                        <bool>; AND A--C
                        <bool>; AND A--D
                        <bool>; AND B--C""")
    assert isinstance(quick, qui.Quick)
    assert len(quick._reactions) == 1
    assert len(quick._contingencies) == 1

def test_quick_multiple_reactions_contingencies():
    # @basti idem
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
    # @basti This is probably not finished :) I'll try to fix this parsing error.
    quick = qui.Quick("A_ppi_B; ! A_[n]--[m]")
    pass