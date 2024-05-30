import CMGDB
import dytop.RoA as RoA


import numpy as np
import matplotlib.pyplot as plt  # module to plot graphs

if __name__ == "__main__":
    # select how many subdivisions will be performed to discretize the state space into boxes
    phase_subdiv = 3

    # setup the bounds for the discretization
    lower_bounds = [-3.0]
    upper_bounds = [3.0]

    # the discretization of the state space will produce a grid

    # create a multivalued combinatorial map


    def f(x):
        return [np.arctan(4*x[0])]

    # Define box map for f


    def F(rect):
        x_min, x_max = CMGDB.BoxMap(f, rect, padding=False)
        d = (rect[1]-rect[0])/4
        return [x_min - d, x_max + d]


    # CMGDB
    model = CMGDB.Model(phase_subdiv, lower_bounds, upper_bounds, F)


    morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)

    roa = RoA.RoA(map_graph, morse_graph)

    roa.save_file()

    fig, ax = roa.PlotRoA()
    plt.show()