import pytest
from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_graph.reaction_graph import ReactionGraph

def test_simple_system():
    quick_system = Quick("A_[b]_ppi+_B_[a]")
    reaction_graph = ReactionGraph(quick_system.rxncon_system)

