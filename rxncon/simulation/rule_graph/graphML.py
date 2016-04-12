import typing as tp
import networkx as nex


G = nex.Graph()
G.add_nodes_from(["a", "b", "c"])
#G.add_node('d', {'domain':"d", 'subdomain': "s", 'residue': 'r'})
G.add_node('e', dict(name="e", id='d', value=30))


edges = [("a", "b"), ("a", "c"), ('b', 'e')]

G.add_edges_from(edges)

print(G.number_of_edges())
print(G.number_of_nodes())
print(G.number_of_selfloops())

gml = nex.generate_gml(G)
nex.write_graphml(G, '/home/thiemese/project/rxncon/graphml/test1.gml')


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
        for onenode in self.graph.nodes(data=True):
            id = onenode[0]
            attr = dict(onenode[1])

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

        for oneedge in self.graph.edges(data=True):
            edge = '<edge source="{}" target="{}">'.format(oneedge[0], oneedge[1])
            for k, v in oneedge[2].items():
                edge += '<att name="{}" value="{}" type="string" />'.format(k, v)
            edge += '</edge>'
            edges.append(edge)
        return "\n".join(edges)

    def xgmmlwriter(self, file: str):

        with open(file, "w") as writehandle:
            writehandle.write(self.to_xgmml())


test_graph = XGMML(G, 'test_graph')
print(test_graph.to_xgmml())
test_graph.xgmmlwriter("/home/thiemese/project/rxncon/graphml/test.xgmml")