import os


class XGMML:
    def __init__(self, graph, graph_name: str):
        self.graph = graph
        self.graph_name = graph_name

    def to_string(self):
        xgmml = [self._header_string(), self._nodes_string(), self._edges_string(), self._footer_string()]
        return "\n".join(xgmml)

    def to_file(self, file_path):
        path, file = os.path.split(file_path)
        if path and os.path.exists(path):
            if not os.path.isfile(file_path):
                self._write_to_file(file_path)
            else:
                raise FileExistsError("{0} exists! remove file and run again".format(file_path))
        elif not path:
            if not os.path.isfile(file):
                self._write_to_file(file_path)
            else:
                print(os.path.dirname(file_path))
                raise FileExistsError("{0} exists! remove file and run again".format(file_path))
        elif path and not os.path.exists(path):
            raise NotADirectoryError("Path {0} does not exists.".format(path))

    def _write_to_file(self, file_path):
        with open(file_path, mode='w') as writehandle:
            writehandle.write(self.to_string())

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
