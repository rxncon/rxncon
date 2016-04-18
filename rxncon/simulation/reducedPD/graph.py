import typing as tp
from enum import unique
from rxncon.util.utils import OrderedEnum
import rxncon.simulation.rule_based.rule_based_model as rbm

@unique
class Directionality(OrderedEnum):
    unidirectional = "->"
    bidirectional = "<->"


class Node:
    pass


class ReactantNode(Node):

    def __init__(self, reactant: rbm.Reactant):
        self.reactant = reactant
#        self.children = []

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{}".format(self.reactant)

    def __eq__(self, other):
        return self.reactant == other.reactant


class ReactionNode(Node):
    def __init__(self, reaction_type: Directionality):
        self.reaction_type = reaction_type

    def __str__(self):
        return "{}".format(self.reaction_type.name)


class Edge:
    def __init__(self, source: Node, target: tp.Union['Node', 'Edge'], type: EdgeType):
        self.source = source
        self.target = target

    def __str__(self):
        return "souce: {0} target: {1}".format(self.source, self.target, self.type.name)


class Graph:
    def __init__(self):
        self.nodes = []
        self.edges = []

    def add_node(self, node: Node):
        self.nodes.append(node)

    def remove_node(self, node: Node):
        self.nodes.remove(node)

    def replace_node(self, node: Node):
        pass

    def add_edge(self, edge: Edge):
        self.edges.append(edge)

    def remove_edge(self, edge: Edge):
        self.edges.remove(edge)

    def neighboring_nodes(self, node: Node):
        neighboring_nodes = [edge.target for edge in self.edges if edge.source == node and isinstance(edge.target, Node)]
        return neighboring_nodes


    def show_nodes(self, node: Node, level = 0):
        if node  in self.nodes:
            neighboring_nodes = g.neighboring_nodes(node)
            print('{0} {1}'.format("\t" * level, node))
            if neighboring_nodes:
                level += 1
                for node in neighboring_nodes:
                    self.show(node, level)
        else:
            raise NameError("Node {} not in graph!".format(node))

if __name__ == '__main__':
    g = Graph()
    g.add_node(ReactantNode('A'))
    g.add_node(ReactantNode('B'))
    g.add_node(ReactantNode('Bp'))
    g.add_node(ReactionNode(Directionality.unidirectional))

    g.add_edge(Edge(ReactantNode('B'), ReactionNode(Directionality.unidirectional)))
    g.add_edge(Edge(ReactantNode('A'), ReactionNode(Directionality.unidirectional)))
    g.add_edge(Edge(ReactionNode(Directionality.unidirectional), ReactantNode('Bp')))

    g.show(Node('B'))
