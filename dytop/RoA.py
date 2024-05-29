# RoA.py  # 2022-20-01
# MIT LICENSE 2020 Ewerton R. Vieira

import pychomp
import CMGDB

import numpy as np
import csv
import os

import matplotlib
import dytop.Poset as Poset
import dytop.Morse_graph_updated as Morse_graph_updated
import dytop.PlotRoA as PlotRoA


class RoA:

    def propagate_adaptive(self, u, adjacencies):
        """propagate a subtree with root u and
        assign the maximal morse node in the subtree for each tile. Inclui adaptive grid"""

        adjacencies_morse_node = set()  # save morse node assigned to each adjacent tile
        for w in adjacencies:
            morse_node = self.dict_tiles.get(w, None)  # get the morse node assigned to tile

            if morse_node == None:  # if there isnt morse node assigned then propagate
                morse_node = self.propagate_adaptive(w, self.map_graph.adjacencies(w))
            adjacencies_morse_node.add(morse_node)  # add the morse node assigned to tile w

        # get one of the maximal morse node in the subtree, it can have more than one,
        # so here we select the first ()
        # tiles that are mapped outside are assign to a fake node with value equal to -1

        adjacencies_morse_node = adjacencies_morse_node - {-1}
        max = list(self.MG.maximal(adjacencies_morse_node))

        if max:
            morse_node = max[0]

            if len(max) == 1:  # assign only for unique since we cant create more relations
                for new_morse_node in adjacencies_morse_node - self.MG_updated.adjacencies(morse_node):
                    self.MG_updated.add_edge(morse_node, new_morse_node)

        else:
            morse_node = -1

        self.dict_tiles[u] = morse_node
        if u in self.S:  # remove it since we dont have to compute again
            self.S.remove(u)
        return morse_node

    def assign_morse_nodes2tiles_adaptive(self):
        """For each tile assign a morse node (creating a region of attraction for a downset)
         Inclui adaptive grid"""

        # It will create a updated Morse graph
        self.MG_updated = pychomp.DirectedAcyclicGraph()
        for vertice in range(self.morse_graph.num_vertices()):
            self.MG_updated.add_vertex(vertice)


        # clone self.tiles_in_morse_sets we dont have to recompute
        self.dict_tiles = dict(self.tiles_in_morse_sets)

        # keep the keys in self.tiles_in_morse_sets we have to recompute
        #TODO: remove keys that correspond to attractors
        self.S = list(self.vertices())

        while self.S:  # loop: assign morse node to a tile and remove it
            v = self.S[0]
            # propagate to determine which morse node should be assign to tile v
            morse_node = self.propagate_adaptive(v, self.map_graph.adjacencies(v))


        # self.MG_updated = self.MG_updated.transitive_reduction
        self.morse_graph_updated = Morse_graph_updated.Morse_graph_updated(self.morse_graph, self.MG_updated)

        return self.dict_tiles

    def propagate(self, u, adjacencies):
        """propagate a subtree with root u and
        assign the maximal morse node in the subtree for each tile"""

        adjacencies_morse_node = set()  # save morse node assigned to each adjacent tile
        for w in adjacencies:
            morse_node = self.dict_tiles.get(w, None)  # get the morse node assigned to tile

            if morse_node == None:  # if there isnt morse node assigned then propagate
                morse_node = self.propagate(w, self.map_graph.adjacencies(w))
            adjacencies_morse_node.add(morse_node)  # add the morse node assigned to tile w

        # get one of the maximal morse node in the subtree, it can have more than one,
        # so here we select the first ()
        # tiles that are mapped outside are assign to a fake node with value equal to -1

        max = list(self.MG.maximal(adjacencies_morse_node - {-1}))
        morse_node = max[0] if max else -1

        self.dict_tiles[u] = morse_node
        if u in self.S:  # remove it since we dont have to compute again
            self.S.remove(u)
        return morse_node

    def assign_morse_nodes2tiles(self):
        """For each tile assign a morse node (creating a region of attraction for a downset)"""

        # clone self.tiles_in_morse_sets we dont have to recompute
        self.dict_tiles = dict(self.tiles_in_morse_sets)

        # remove keys in self.tiles_in_morse_sets we dont have to recompute
        self.S = list(self.vertices() - set(self.tiles_in_morse_sets.keys()))

        while self.S:  # loop: assign morse node to a tile and remove it
            v = self.S[0]
            # propagate to determine which morse node should be assign to tile v
            morse_node = self.propagate(v, self.map_graph.adjacencies(v))

        return self.dict_tiles

    def __init__(self, map_graph, morse_graph, adaptive=False):
        """
        Region of Attraction class
        Assign cells in the phase space that are mapped to a unique Morse Node
        (Regions that are uniquely mapped to Morse Sets).
        Equivalent to an order retraction onto the Morse tiles by mapping
        to unique successor.
        Input: adaptive = True: Adaptive grid and update Morse graph
        """

        self.dir_path = os.path.abspath(os.getcwd()) + "/output/"

        self.morse_graph = morse_graph
        self.map_graph = map_graph

        # Get number of vertices
        self.num_verts = map_graph.num_vertices()

        self.vertices_ = {a for a in range(self.num_verts)}

        self.tiles_in_morse_sets = {}

        # create a dict: tile in morse set -> morse node
        # it is need to create the condensation graph
        for i in range(self.morse_graph.num_vertices()):
            for j in self.morse_graph.morse_set(i):

                self.tiles_in_morse_sets[j] = i

        cyclic_Morse_graph = False
        MG = .DirectedAcyclicGraph()  # building Morse Graph Poset
        MG.add_vertex(0)
        for u in range(morse_graph.num_vertices()):
            for v in morse_graph.adjacencies(u):
                MG.add_edge(u, v)
                if u in MG.adjacencies(v):  # prevent to compute with cyclic MG
                    cyclic_Morse_graph = True

        if cyclic_Morse_graph:
            print("\033[91m Morse Graph input is cyclic, wrong input \033[0m")
            self.dict_tiles = False

        else:
            self.MG = Poset.Poset(MG)


            if adaptive:
                self.assign_morse_nodes2tiles_adaptive()
            else:
                self.assign_morse_nodes2tiles()

        # print(
        #     f"memory size of: tiles_in_morse_sets={asizeof(self.tiles_in_morse_sets)}\n",
        #     f"memory size of: dict_tiles={asizeof(self.dict_tiles)}\n",
        #     f"MG={asizeof(self.MG)}")

        # newcolors
        viridis = matplotlib.cm.get_cmap('viridis', 256)
        newcolors = viridis(np.linspace(0, 1, 256))
        orange = np.array([253/256, 174/256, 97/256, 1])
        yellowish = np.array([233/256, 204/256, 50/256, 1])
        newcolors[109:146, :] = orange
        newcolors[219:, :] = yellowish
        self.newcmp = matplotlib.colors.ListedColormap(newcolors)
        # self.newcmp=matplotlib.cm.brg  # old default option

    def vertices(self):
        """
        Return the set of elements in the poset
        """
        return self.vertices_

    def box_center(self, rect):
        dim = len(rect) // 2
        return [rect[i] + (rect[dim + i] - rect[i])/2 for i in range(dim)]

    def save_file(self, name=""):
        rect = self.morse_graph.phase_space_box(0)
        dim = int(len(rect) // 2)
        size_box = [rect[dim + i] - rect[i] for i in range(dim)]

        # getting the bounds
        lower_bounds = rect[0:dim]
        rect = self.morse_graph.phase_space_box(self.map_graph.num_vertices()-1)
        upper_bounds = rect[dim::]

        name = self.dir_path + name + "_RoA_" + ".csv"
        with open(name, "w") as file:
            f = csv.writer(file)
            f.writerow(["Box size", "Lower bounds", "Upper bounds"])
            f.writerow(size_box+lower_bounds+upper_bounds)
            f.writerow(["Tile", "Morse_node", "Box"])

            if self.dict_tiles:  # it might be empty
                # tiles in roa
                for tile_in_roa in set(self.dict_tiles.items()) - set(self.tiles_in_morse_sets.items()):
                    tile_in_roa = list(tile_in_roa) + \
                        [a for a in self.morse_graph.phase_space_box(tile_in_roa[0])]
                    f.writerow(tile_in_roa)
                # tiles in morse sets
                f.writerow(["Tile_in_Morse_set", "Morse_node", "Box"])
                for tile_in_morse_set in self.tiles_in_morse_sets.items():
                    tile_in_morse_set = list(tile_in_morse_set) + \
                        [a for a in self.morse_graph.phase_space_box(tile_in_morse_set[0])]
                    f.writerow(tile_in_morse_set)

    def Morse_sets_vol(self):
        """Compute a dict that gives the volume of the regions of attraction"""

        d_vol = dict()
        tiles_and_morse_nodes = list(self.dict_tiles.items())
        for tile_and_morse in tiles_and_morse_nodes:
            i, j = tile_and_morse  # i is the tile and j is the associated morse node
            box = self.morse_graph.phase_space_box(i)

            size = len(box)
            half = int(size / 2)

            volume_cube = 1
            for k in range(half):
                volume_cube *= float(box[half + k]) - float(box[k])

            d_vol[j] = d_vol.get(j, 0) + volume_cube

        return d_vol

    def PlotRoA(self, selection=[], fig_w=8, fig_h=8, xlim=None, ylim=None,
                  cmap=matplotlib.cm.get_cmap('viridis', 256), name_plot='', from_file=None, plot_point=False, section=None):

        if not os.path.exists(self.dir_path):
            os.makedirs(self.dir_path)
        self.save_file(name="temp")

        rect = self.morse_graph.phase_space_box(0)
        dim = int(len(rect) // 2)

        # getting the bounds
        lower_bounds = rect[0:dim]
        upper_bounds = rect[dim::]

        for i in range(self.num_verts):
            box = self.morse_graph.phase_space_box(i)
            for index, j in enumerate(box[0:dim]):
                if lower_bounds[index] > j:
                    lower_bounds[index] = j
            for index, j in enumerate(box[dim::]):
                if upper_bounds[index] < j:
                    upper_bounds[index] = j

        fig, ax = PlotRoA.PlotRoA(lower_bounds, upper_bounds, selection=selection, fig_w=fig_w, fig_h=fig_h, xlim=xlim,
                            ylim=ylim, cmap=cmap, name_plot=name_plot, from_file="temp", plot_point=plot_point, section=section)

        os.remove(self.dir_path + "temp_RoA_.csv")
        return fig, ax

