# CMGDB_util.py  # 2021-10-26
# MIT LICENSE 2020 Ewerton R. Vieira


import matplotlib
import numpy as np
import csv
import os
import CMGDB
from datetime import datetime


def base10toN(num,n):
    """N<10"""
    new_num_string=''
    current=num
    while current!=0:
        remainder=current%n
        remainder_string=str(remainder)
        new_num_string=remainder_string+new_num_string
        current=current//n
    return new_num_string





# rect = [0,0,2,2]
# subdivisions = 2
#
# print(len(FacePoints(rect, subdivisions)))
# print(FacePoints(rect, subdivisions))


class CMGDB_util:

    def __init__(self):
        self.dir_path = os.path.abspath(os.getcwd()) + "/output/"

    def Local_Lipschitz(self, X, Y):
        """Computes an estimation of the local Lipschitz constant of a given set of points"""
        dim = len(X[0])
        slope_max = [0 for i in range(dim)]
        while len(X) > 1:
            x0 = X.pop()
            x0 = np.array(x0)
            y0 = Y.pop()
            y0 = np.array(y0)
            for k, x1 in enumerate(X):
                x1 = np.array(x1)
                y1 = np.array(Y[k])
                for i in range(dim):
                    dy = np.linalg.norm(y1[i] - y0[i])
                    if np.linalg.norm(y1[i] - y0[i]) > 3:
                        continue


                    slope = dy / np.linalg.norm(x1 - x0)
                    if slope > slope_max[i]:
                        slope_max[i] = slope
        return slope_max

    def FacePoints(self, rect, subdivisions):
        """return the mid and the corners points of all faces of a rectangle"""
        dim = len(rect) // 2
        low = rect[0:dim]
        up = rect[dim::]
        box_size = [(rect[i+dim] - rect[i])/subdivisions for i in range(dim)]
        diagonal = [[rect[i] + j*box_size[i] for i in range(dim)] for j in range(1, subdivisions)]


        all = []

        for face_index in range(2**dim -1):  # all faces, remember that 2**dim is inside the cube
            face = [face_index][0]  # copy to avoid changing the indexing
            position = 0  # postion of the bits of face

            face_points = []
            while face or position < dim:
                face_points_temp = face_points.copy()   # copy to avoid changing the indexing
                face_points = []
                if not face_points_temp:
                    face_points_temp = [[]]

                if face & 1:  # if 1 than subdivision
                    for point in face_points_temp:
                        for diag_element in diagonal:
                            face_points.append(point + [diag_element[position]])

                else:  # add low or up elements of the bounds
                    for point in face_points_temp:
                        face_points.append(point + [low[position]])
                        face_points.append(point + [up[position]])

                face = face >> 1  # next bit
                position += 1  # next bit position

            all += face_points

        return all

    def FacePoints_old(self, rect):
        """return the mid and the corners points of all faces of a rectangle"""
        dim = len(rect) // 2
        face = rect + [(rect[i+dim] + rect[i])/2 for i in range(dim)]
        all = []
        for i in range(3**dim):
            ii = base10toN(i,3)
            ii = ii[::-1] #invert
            ii += '0'*(dim - len(ii)) #add 0
            all.append([face[j + int(ii[j])*dim] for j in range(dim)])
        return all

    def run_CMGDB(self, subdiv_min, subdiv_max, lower_bounds, upper_bounds, phase_periodic, F, base_name, subdiv_init=6, subdiv_limit=10000, cmap=matplotlib.cm.get_cmap('viridis', 256)):
        # Define the parameters for CMGDB

        model = CMGDB.Model(subdiv_min, subdiv_max, subdiv_init, subdiv_limit,
                            lower_bounds, upper_bounds, phase_periodic, F)

        startTime = datetime.now()

        morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)

        print(datetime.now() - startTime)

        # Save Morse graph
        
        MG = self.dir_path + base_name
        morse_fname = self.dir_path + base_name + ".csv"

        CMGDB.PlotMorseGraph(morse_graph, cmap).format = 'png'
        CMGDB.PlotMorseGraph(morse_graph, cmap).render(MG)

        # Save file
        morse_nodes = range(morse_graph.num_vertices())
        morse_sets = [box + [node]
                      for node in morse_nodes for box in morse_graph.morse_set_boxes(node)]
        np.savetxt(morse_fname, np.array(morse_sets), delimiter=',')

        return morse_graph, map_graph

    # def run_CMGDB(phase_subdiv, lower_bounds, upper_bounds, phase_periodic, F, base_name):
    #     # Define the parameters for CMGDB
    #
    #     model = CMGDB.Model(phase_subdiv, lower_bounds, upper_bounds, phase_periodic, F)
    #
    #     startTime = datetime.now()
    #
    #     morse_graph, map_graph = CMGDB.ComputeMorseGraph(model)
    #
    #     print(datetime.now() - startTime)
    #
    #     # Save Morse graph
    #     self.dir_path = os.path.abspath(os.getcwd()) + "/output/"
    #     MG = self.dir_path + base_name
    #     morse_fname = self.dir_path + base_name + ".csv"
    #
    #     CMGDB.PlotMorseGraph(morse_graph).format = 'png'
    #     CMGDB.PlotMorseGraph(morse_graph).render(MG)
    #
    #     # Save file
    #     morse_nodes = range(morse_graph.num_vertices())
    #     morse_sets = [box + [node] for node in morse_nodes for box in morse_graph.morse_set_boxes(node)]
    #     np.savetxt(morse_fname, np.array(morse_sets), delimiter=',')
    #
    #     return morse_graph, map_graph

    def F_K(self, f, rect, K):
        """Input: function f, rectangle rect, and the Lipschit constant in vector form K
        Output: Image of rectangle by f, also taking in account the expansion given by K"""
        half = len(rect) // 2
        im_center = f([rect[i] + (rect[half + i] - rect[i])/2 for i in range(half)])
        list1 = []
        list2 = []
        for i in range(half):  # image of the center of the rect +or- lenght * K
            list1.append(im_center[i] - (rect[half + i] - rect[i]) * K[i] / 2)
            list2.append(im_center[i] + (rect[half + i] - rect[i]) * K[i] / 2)
        return list1 + list2

    def Morse_sets_vol(self, name_file):
        """Compute the volume of a Morse set"""
        with open(name_file, 'r') as f:
            lines = csv.reader(f, delimiter=',')
            d_vol = dict()
            for row in lines:
                size = len(row) - 1
                half = int(size / 2)
                volume_cube = 1
                for i in range(half):
                    volume_cube *= float(row[half + i]) - float(row[i])
                if row[size] in d_vol.keys():
                    d_vol[row[size]] = d_vol[row[size]] + volume_cube
                else:
                    d_vol[row[size]] = volume_cube
        return d_vol

    def sample_points(self, lower_bounds, upper_bounds, num_pts):
        # Sample num_pts in dimension dim, where each
        # component of the sampled points are in the
        # ranges given by lower_bounds and upper_bounds
        dim = len(lower_bounds)
        X = np.random.uniform(lower_bounds, upper_bounds, size=(num_pts, dim))
        return X

    def BoxMapK(self, f, rect, K):
        dim = int(len(rect) / 2)
        X = CMGDB.CornerPoints(rect)
        # Evaluate f at point in X
        Y = [f(x) for x in X]
        # Get lower and upper bounds of Y
        Y_l_bounds = [min([y[d] for y in Y]) - K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
        Y_u_bounds = [max([y[d] for y in Y]) + K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
        return Y_l_bounds + Y_u_bounds


    def BoxMapK_valid(self, f, rect, K, valid_list, point2cell):
        """Box map for valid cells based on a list using point2cell indexing"""
        id = point2cell(CMGDB.CenterPoint(rect)[0]) # center to avoid boundary points
        dim = len(rect) // 2
        if valid_list[id]:
            X = self.FacePoints(rect, 2)
            # X = CMGDB.CornerPoints(rect)
            Y = [f(x) for x in X]
            # Get lower and upper bounds of Y
            Y_l_bounds = [min([y[d] for y in Y]) - K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
            Y_u_bounds = [max([y[d] for y in Y]) + K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
            return Y_l_bounds + Y_u_bounds
        else:
            return [30000]*2*dim

    def BoxMapK_AE(self, f, rect, K, valid_state):
        """Box map for autoencoder"""
        dim = int(len(rect) / 2)
        X = CMGDB.CornerPoints(rect)
        # Evaluate f at point in X
        # print(X)
        Y = [f(x) for x in X if valid_state(x)]
        # print(Y)

        # print(Y)

        # if not Y: return [0]*2*dim
        # if len(Y)<4: return [3]*2*dim#[-.5, .25, -.5+0.01, .25+0.01]#[0]*2*dim
        if len(Y)<4 or np.linalg.norm(np.array(X[0]))<0 : return [3]*2*dim#[-.5, .25, -.5+0.01, .25+0.01]#[0]*2*dim

        # Get lower and upper bounds of Y
        Y_l_bounds = [min([y[d] for y in Y]) - K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
        Y_u_bounds = [max([y[d] for y in Y]) + K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
        return Y_l_bounds + Y_u_bounds


    def Box_GP_K(self, learned_f, rect, K, n=-3):
        """learned_f with predicted mean and standard deviation
        K Lipschit constant"""
        dim = int(len(rect) / 2)
        X = CMGDB.CornerPoints(rect)
        # print(X)
        # Evaluate f at point in X

        Y, S = learned_f(X[0])

        X = X[1::]
        for x_ in X:
            y_, s_ = learned_f(x_)
            Y = np.concatenate((Y, y_))
            S = np.concatenate((S, s_))

        # print(f"{Y}\n {S} \n {S_}")
        Y_max = Y + S*(2**n)
        Y_min = Y - S*(2**n)

        # print(f"{Y_min}\n {Y_max}")

        # print(f"{np.amin(Y_min[:,0])}\n {np.amax(Y_max[:,1])}")

        # Get lower and upper bounds of Y
        Y_l_bounds = [np.amin(Y_min[:, d]) - K*(rect[d + dim] - rect[d]) for d in range(dim)]
        Y_u_bounds = [np.amax(Y_max[:, d]) + K*(rect[d + dim] - rect[d]) for d in range(dim)]
        return Y_l_bounds + Y_u_bounds

    def F_GP_K(self, learned_f, rect, K, n=-3):
        """learned_f with predicted mean and standard deviation
        K Lipschit constant"""
        dim = int(len(rect) / 2)
        X = CMGDB.CenterPoint(rect)
        # print(X)
        # Evaluate f at point in X
        Y, S = learned_f(X)

        # print(f"{Y}\n {S} \n {S_}")
        Y_max = Y + S*(2**n)
        Y_min = Y - S*(2**n)

        # print(f"{Y_min}\n {Y_max}")

        # print(f"{np.amin(Y_min[:,0])}\n {np.amax(Y_max[:,1])}")

        # Get lower and upper bounds of Y
        Y_l_bounds = [np.amin(Y_min[:, d]) - K*(rect[d + dim] - rect[d]) for d in range(dim)]
        Y_u_bounds = [np.amax(Y_max[:, d]) + K*(rect[d + dim] - rect[d]) for d in range(dim)]
        return Y_l_bounds + Y_u_bounds

    def Box_J(self, f, J, rect):
        """f: function, J: Jacobian matrix, rect: rectangle
        Given a rect return the smallest rectangle that contains the image of the
        J \cdot rect"""
        dim = int(len(rect) / 2)
        x = rect[0:dim]
        y = f(x)
        Jac, _ = J(np.array(x).reshape(-1, dim))
        # next, add the sum of the columns
        Jac = np.concatenate((Jac, Jac.sum(axis=0).reshape(1, dim)), axis=0)

        Y_l_bounds = []
        Y_u_bounds = []
        for d in range(dim):
            Y_l_bounds.append(
                np.amin(Jac[:, d]) * (rect[d + dim] - rect[d]) + y[d]
            )

            Y_u_bounds.append(
                np.amax(Jac[:, d]) * (rect[d + dim] - rect[d]) + y[d]
            )

        return Y_l_bounds + Y_u_bounds

    def F_J(self, learned_f, J, rect, lower_bounds, n=-3, weak_multivalued_map=1):
        """f: function, J: Jacobian matrix, rect: rectangle
        Given a rect return the smallest rectangle that contains the image of the
        J \cdot rect
        weak_multivalued_map: the number to compute the subset S of all cells"""
        dim = int(len(rect) / 2)
        X = rect[0:dim]
        Y, S = learned_f(X)
        Y_max = Y + S*(2**n)
        Y_min = Y - S*(2**n)

        size_of_box = [rect[d+dim] - rect[d] for d in range(dim)]

        if weak_multivalued_map:
            coordinate = [int(
                np.rint((rect[d] - lower_bounds[d]) / size_of_box[d])
            ) % (1+weak_multivalued_map) for d in range(dim)
            ]
        else:
            coordinate = True

        # print(coordinate)
        Y_l_bounds = []
        Y_u_bounds = []

        # print(f"X {X} \n Y {Y} \n Y_min {Y_min}\n Y_max {Y_max} \n size_of_box {size_of_box}")

        if not all(coordinate):  # Compute J if all coordinate = 0 module (weak_multivalued_map+1)

            print(coordinate)
            Jac, _ = J(np.array(X).reshape(-1, dim))

            # print(f"J {Jac}")
            for d in range(dim):
                J_d_norm = np.linalg.norm(Jac[d, :])
                Y_l_bounds.append(Y_min[:, d] - J_d_norm * size_of_box[d])
                Y_u_bounds.append(Y_max[:, d] + J_d_norm * size_of_box[d])

        else:

            for d in range(dim):
                Y_l_bounds.append(Y_min[:, d])
                Y_u_bounds.append(Y_max[:, d])

        return Y_l_bounds + Y_u_bounds

    def Box_ptwise(self, learned_f, rect, n=-3):
        """learned_f with predicted mean applied to the corner points
        and standard deviation applied to the center point"""

        dim = int(len(rect) / 2)
        X = CMGDB.CenterPoint(rect) + CMGDB.CornerPoints(rect)
        # Evaluate f at point in X

        Y, S = learned_f(X[0])

        Y = Y - S*(2**n)
        Y = np.concatenate((Y, Y + 2 * S * (2 ** n)))

        X = X[1::]
        for x_ in X:
            y_, _ = learned_f(x_)
            Y = np.concatenate((Y, y_))

        # Get lower and upper bounds of Y
        Y_l_bounds = [np.amin(Y[:, d]) for d in range(dim)]
        Y_u_bounds = [np.amax(Y[:, d]) for d in range(dim)]
        return Y_l_bounds + Y_u_bounds

    def Box_noisy_K(self, f, rect, K, noise):
        """Box map with noise = [x_epsilon, y_epsilon, z_epsilon, f_epsilon]"""

        noise_x = noise[0:-1]
        noise_f = noise[-1]

        dim = int(len(rect) / 2)
        X = CMGDB.CornerPoints(rect)
        # X = self.FacePoints(rect, 2)
        # Evaluate f at point in X
        Y = [f(x) for x in X]

        # K = 2 * self.Local_Lipschitz(X, Y)

        # Get lower and upper bounds of Y
        Y_l_bounds = [min([y[d] for y in Y]) - K[d] * (rect[d + dim] - rect[d] + noise_x[d]) - noise_f for d in range(dim)]
        Y_u_bounds = [max([y[d] for y in Y]) + K[d] * (rect[d + dim] - rect[d] + noise_x[d]) + noise_f for d in range(dim)]
        f_rect = Y_l_bounds + Y_u_bounds
        return f_rect

    def Box_noisy_K_wa(self, f, rect, K, noise):
        """Box map with noise = [x_epsilon, y_epsilon, z_epsilon, f_epsilon]"""

        noise_x = noise[0:-1]
        noise_f = noise[-1]

        dim = int(len(rect) / 2)
        X = CMGDB.CornerPoints(rect)
        # X = self.FacePoints(rect, 2)
        # Evaluate f at point in X
        Y = np.array([f(x) for x in X])

        # K = 2 * self.Local_Lipschitz(X, Y)
        
        Amin = np.amin(Y, axis=0)
        Amax = np.amax(Y, axis=0)
        Y_l_bounds = [Amin[d] - K[d] * (rect[d + dim] - rect[d] + noise_x[d]) - noise_f for d in range(dim)]
        Y_u_bounds = [Amax[d] + K[d] * (rect[d + dim] - rect[d] + noise_x[d]) + noise_f for d in range(dim)]
        if Amax[0] - Amin[0]>3.14:
            # print("DIFF", Amax[0] - Amin[0])
            Y_l_bounds[0] = Y_u_bounds[0]
            Y_u_bounds[0] = 3.14158
        f_rect = Y_l_bounds + Y_u_bounds
        return f_rect

    def F_data(self, rect, id2image, point2cell, K):
        dim = len(rect) // 2
        id_of_rect = point2cell(CMGDB.CenterPoint(rect)[0]) # center to avoid boundary points
        Y = id2image[id_of_rect]
        if Y:
            # Y is a list of a array
            # Get lower and upper bounds of Y
            Amin = np.amin(Y, axis=0)
            Amax = np.amax(Y, axis=0)
            Y_l_bounds = [Amin[d] - K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
            Y_u_bounds = [Amax[d] + K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
            f_rect = Y_l_bounds + Y_u_bounds
            return f_rect
        else:
            return [30000]*2*dim

    def F_data_wa(self, rect, id2image, point2cell, K):
        dim = len(rect) // 2
        center = CMGDB.CenterPoint(rect)[0]
        # box_removed = [-0.1, 0.07, 0.1, 0.2]
        # if all([box_removed[d] < center[d] < box_removed[d+dim] for d in range(dim)]):
        #     return [30000]*2*dim
        # box_removed = [-1.5, -0.1, -0.1, 1.6]
        # if all([box_removed[d] < center[d] < box_removed[d+dim] for d in range(dim)]):
        #     return [30000]*2*dim
        # box_removed = [0.1, -0.1, 1.5, 1.6]
        # if all([box_removed[d] < center[d] < box_removed[d+dim] for d in range(dim)]):
        #     return [30000]*2*dim           
        id_of_rect = point2cell(center) # center to avoid boundary points
        Y = id2image[id_of_rect]
        

        if Y:
            # Y is a list of a array
            # Get lower and upper bounds of Y
            Amin = np.amin(Y, axis=0)
            Amax = np.amax(Y, axis=0)
            Y_l_bounds = [Amin[d] - K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
            Y_u_bounds = [Amax[d] + K[d]*(rect[d + dim] - rect[d]) for d in range(dim)]
            if Amax[0] - Amin[0]>3.14:
                # print("DIFF", Amax[0] - Amin[0])
                Y_l_bounds[0] = Y_u_bounds[0]
                Y_u_bounds[0] = 3.14159


            f_rect = Y_l_bounds + Y_u_bounds
            return f_rect
        else:
            return [30000]*2*dim