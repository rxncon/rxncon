import pytest
import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_graph.regulatory_graph as reg
import rxncon.simulation.rule_graph.graphML as gml
def test_regulatory_graph():
    quick_system = qui.Quick('''A_ppi_B; ! A-{P}
                                C_p+_A''')

    reg_system = reg.RegulatoryGraph(quick_system.rxncon_system)
    graph = reg_system.to_graph()
    xgmml_system = gml.XGMML(graph, 'test_graph')
    print(xgmml_system.to_xgmml())
    xgmml_system.xgmmlwriter("/home/thiemese/project/rxncon/graphml/test.xgmml")

