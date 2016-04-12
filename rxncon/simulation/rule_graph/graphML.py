import typing as tp
import networkx as nex
import rxncon.core.rxncon_system as rxs
import rxncon.core.specification as spec
# for cytoscape export:
# reaction graph of rxncon: visalisation of specifications and there relationships
# for each spec one node of certain size
# for each domain/subdomain/residue one node of certain size with edge length 0 to it's main node,
# the length of a edge should be defined by it's weight since there is no algorithm considering the length, which makes
# the length of an edge dependent on the node position only

def _create_node(specification)
def create_detailed_reaction_graph(rxncon_system: rxs.RxnConSystem):
    G = nex.Graph()
    for reaction in rxncon_system.reactions:
        _create_node(reaction.object)
        _create_node(reaction.subject)

# G.add_nodes_from(["a", "b", "c"])
# #G.add_node('d', {'domain':"d", 'subdomain': "s", 'residue': 'r'})
# #G.add_node('e', dict(name="e", size=40))
# #G.add_node('f', dict(name="f", size=20))
#
# G.add_node('e', dict(name='e', size=20))
# G.add_node('f')
#
#
# G.add_edge('b', 'e', weight=10, length=10)
# G.add_edge("a", "c", weight=10, length=10)
# G.add_edge("a", "b", weight=10, length=10)
# G.add_edge('e', 'f', weight=0.0, length=0.0)


class XGMML:
    def __init__(self, graph, graph_name: str):
        self.graph = graph
        self.graph_name = graph_name

    def to_xgmml(self):
        xgmml = [self._header_string(), self._nodes_string(), self._edges_string(), self._footer_string()]
        return "\n".join(xgmml)

    def _header_string(self):
        return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
    <att name="selected" value="1" type="boolean" />
    <att name="name" value="{0}" type="string"/>
    <att name="shared name" value="{0}" type="string"/>
    """.format(self.graph_name)

    def _footer_string(self):
        return '</graph>'
    def _nodes_string(self):
        nodes = []
        for graph_node in self.graph.nodes(data=True):
            id = graph_node[0]
            attr = dict(graph_node[1])

            if 'label' in attr:
                label = attr['label']
                del attr['label']
            else:
                label = id

            node = '<node id="{id}" label="{label}">'.format(id=id, label=label)
            for k, v in attr.items():
                node += '<att name="{}" value="{}" />'.format(k, v)
            node += '</node>'
            nodes.append(node)
        return "\n".join(nodes)

    def _edges_string(self):
        edges = []

        for graph_edge in self.graph.edges(data=True):
            edge = '<edge source="{}" target="{}">'.format(graph_edge[0], graph_edge[1])
            for k, v in graph_edge[2].items():
                edge += '<att name="{}" value="{}"/>'.format(k, v)
            edge += '</edge>'
            edges.append(edge)
        return "\n".join(edges)

    def xgmmlwriter(self, file: str):

        with open(file, "w") as writehandle:
            writehandle.write(self.to_xgmml())


test_graph = XGMML(G, 'test_graph')
print(test_graph.to_xgmml())
test_graph.xgmmlwriter("/home/thiemese/project/rxncon/graphml/test.xgmml")