# MultivaluedMap.py  # 2021-20-10
# MIT LICENSE 2020 Ewerton R. Vieira
import pychomp
import graphviz
import dytop.Poset as Poset


def Dict_inverse(D):
    """Return the inverse of a dictionary D"""
    inv_map = {}
    for k, v in D.items():
        for i in v:
            inv_map[i] = inv_map.get(i, []) + [k]
    return inv_map


class MultivaluedMap:

    def __init__(self, map_graph, morse_graph):
        """
        Adding features to map_graph to become a mvm_graph
        """

        # self.dir_path = os.path.abspath(os.getcwd()) + "/output/"

        self.morse_graph = morse_graph

        # Get number of vertices
        num_verts = map_graph.num_vertices()

        self.vertices_ = {a for a in range(num_verts)}

        self.map_graph = map_graph

        # self.mapping: map_graph -> condensation_graph
        condensation_graph, self.mapping = pychomp.CondensationGraph(
            self.vertices_, map_graph.adjacencies)

        self.CG = condensation_graph
        # self.CG.graphviz = self.condensation_graph_with_labels

        self.mapping_inv = {}  # inverse of mapping condensation_graph -> map_graph
        for key, value in self.mapping.items():
            self.mapping_inv.setdefault(value, []).append(key)

        # condensation_graph with poset structure
        self.condensation_graph = Poset.Poset(condensation_graph, save_memory=True)

        MG = pychomp.DirectedAcyclicGraph()  # building Morse Graph Poset
        MG.add_vertex(0)
        for u in range(morse_graph.num_vertices()):
            for v in morse_graph.adjacencies(u):
                MG.add_edge(u, v)
        self.MG = Poset.Poset(MG)

    def vertices(self):
        """
        Return the set of elements in the poset
        """
        return self.vertices_

    def upset(self, v):
        """
        Return the set { u : u >= v }
        """
        U = []
        for a in self.condensation_graph.ancestors(self.mapping[v]).union({self.mapping[v]}):
            U = U + self.mapping_inv[a]
        return set(U)

    def mvm_graph(self):
        mvm = pychomp.DirectedAcyclicGraph()  # building domain graph
        mvm.add_vertex(0)
        for v in self.vertices():
            for u in self.map_graph.adjacencies(v):
                mvm.add_edge(v, u)
        return mvm

    def build_viz(self, graph, shape='circle'):
        """ Return a graphviz string describing the graph and its labels """
        gv = 'digraph {\n'

        gv += f'node [fontsize=12, shape={shape}]'

        for v in graph.vertices():
            gv += f"{v}[label={v}];\n"

        for v in graph.vertices():
            for u in graph.adjacencies(v):
                gv += f"{v}->{u};\n"
        return gv + '}\n'

    def condensation_graph_with_labels(self):
        """ Add labels to the condensation_graph"""
        gv = 'digraph {\n'

        gv += f'node [fontsize=12, shape=square]\n'

        for v in self.CG.vertices():
            inv_v = self.mapping_inv[v]
            if len(inv_v) == 1 and inv_v[0] in self.mvm_graph().adjacencies(inv_v[0]):
                gv += f'{v}[label=\"<{inv_v[0]}>\"];\n'
            else:
                temp = ''.join(f'{e} ' for e in self.mapping_inv[v])
                gv += f'{v}[label=\"{temp}\"];\n'

        for v in self.CG.vertices():
            for u in self.CG.adjacencies(v):
                gv += f"{v}->{u};\n"
        gv += "}\n"

        return gv
        # return graphviz.Source(gv)

    def _repr_svg_(self):
        """
        Return svg representation for visual display
        """
        # return graphviz.Source(self.mvm_graph().graphviz())._repr_svg_()

        # f = graphviz.Source(self.graphviz(), filename=self.filename, format='pdf')
        # f.render(cleanup=True);
        # return graphviz.Source(str(self.build_viz(self.mvm_graph())))
        return graphviz.Source(self.build_viz(self.mvm_graph()))
    
    def vizualize_CG(self):
        return graphviz.Source(self.condensation_graph_with_labels())

    def vizualize(self):
        return graphviz.Source(self.build_viz(self.mvm_graph()))
