# Morse_graph_updated.py  # 2023-01-22
# MIT LICENSE 2020 Ewerton R. Vieira

class Morse_graph_updated:
    def __init__(self, morse_graph, MG_updated):

        self.morse_graph = morse_graph

        self.MG_ = Poset.Poset(MG_updated)

        self.MG_updated_ = pychomp.DirectedAcyclicGraph()  # building Morse Graph Poset
        self.MG_updated_.add_vertex(0)

        for u in self.morse_graph.vertices():
            for v in self.MG_.parents(u):
                self.MG_updated_.add_edge(v, u)  # inverse since RoA build MG_updated inverted

        # for (u,v) in self.morse_graph.edges():
        #     if (u,v) in MG_updated.edges():
        #         self.MG_updated_.add_edge(u, v)

    def adjacencies(self, i):
        return self.MG_updated_.adjacencies(i)

    def adjacencies_unreduced(self, i):
        return self.MG_updated_.adjacencies_unreduced(i)

    def annotations(self, i):
        return self.morse_graph.annotations(i)

    def edges(self):
        return self.MG_updated_.edges()

    def edges_unreduced(self):
        return self.MG_updated_.edges_unreduced()

    def morse_set(self, i):
        return self.morse_graph.morse_set(i)

    def morse_set_boxes(self, i):
        return self.morse_graph.morse_set_boxes(i)

    def num_vertices(self):
        return self.morse_graph.num_vertices()

    def phase_space_box(self, i):
        return self.morse_graph.phase_space_box(i)

    def vertices(self):
        return self.morse_graph.vertices()