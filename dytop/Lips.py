# Lips.py  # 2022-11-22
# MIT LICENSE 2020 Ewerton R. Vieira

import scipy.optimize as optimize
import numpy as np
from datetime import datetime


# TODO imcorporate
"""def sampled_Lipschitz(lower_bounds, upper_bounds, N, f, base_name=""):
    x = sample_points(lower_bounds, upper_bounds, N)
    K = 0
    for i in range(len(x)):
        x0 = x[i]
        for j in range(i + 1, len(x)):
            x1 = x[j]
            y0 = np.array(f(x0))
            y1 = np.array(f(x1))
            K_temp = np.linalg.norm(y1 - y0)/np.linalg.norm(x1 - x0)
            if K_temp > K:
                K = K_temp

    dir_path = os.path.abspath(os.getcwd()) + "/output/"
    Lipschitz_file = dir_path + base_name + "_sampled_Lipschitz.txt"
    np.savetxt(Lipschitz_file, [K])

    return [K for a in range(len(lower_bounds))]"""


class Lips:

    def dy_dx(self, x0, x1):
        y0 = np.array(self.g(x0))
        y1 = np.array(self.g(x1))
        x0 = np.array(x0)
        x1 = np.array(x1)
        return np.linalg.norm(y1 - y0), np.linalg.norm(x1 - x0)

    def slope(self, func, x0, x1):
        y0 = np.array(func(x0))
        y1 = np.array(func(x1))
        x0 = np.array(x0)
        x1 = np.array(x1)
        return - np.linalg.norm(y1 - y0)/np.linalg.norm(x1 - x0)

    def f(self, x):
        size = len(x)
        x0 = x[0: size//2]
        x1 = x[size//2: size]
        return self.slope(self.g, x0, x1)

    def Powell(self, x_initial):

        res = optimize.minimize(self.f, x_initial, method='Powell',
                                bounds=self.bnds, options={'xtol': self.error, 'ftol': self.error, 'maxiter': self.N, 'maxfev': self.N})
        # print(res)
        return [res.x, res.fun]

    def acceptance_probability(self, E, E_new, T):
        # if E_new < E:
        #     return 1
        return np.exp(-(E_new-E)/T)

    def random_neighbour(self, x):
        for i in range(len(x)):
            x[i] += np.random.uniform(-1, 1)
            if x[i] < self.bnds[i][0]:
                x[i] = self.bnds[i][0]
            elif x[i] > self.bnds[i][1]:
                x[i] = self.bnds[i][1]
        return x

    def random_x(self):
        return [(self.bnds[i][1] - self.bnds[i][0])*np.random.random() + self.bnds[i][0] for i in range(len(self.bnds))]

    def simulated_annealing(self, steps):
        x = self.random_x()
        X, E = self.Powell(x)
        # print("x=%.2f, fmin=%.2f" % (x, E))
        for k in range(steps):
            T = (1 - k/steps)
            x = self.random_neighbour(x)
            E_new = self.Powell(x)[1]
            X_new = self.Powell(x)[0]
            P = self.acceptance_probability(E, E_new, T)
            if P > np.random.random():
                E = E_new
                X = X_new
            # print("x=%.4f, fmin=%.4f, Prob.=%.4f" % (x, E, P))
        return abs(E), X

    def find_and_save_Lipschitz(self, n):

        startTime = datetime.now()
        E, X = self.simulated_annealing(n)
        print(f"Time to find Lipschitz constant = {datetime.now() - startTime}")

        x_0 = [X[a] for a in range(len(self.bnds)//2)]
        x_1 = [X[a] for a in range(len(self.bnds)//2, len(self.bnds))]
        # print(x_0)
        # print(x_1)
        info = f"Lipschit constant={E}, \nx_0={x_0}, x_1={x_1} \nf(x_0)={self.g(x_0)}, f(x_1)={self.g(x_1)}, \ndy and dx={self.dy_dx(x_0, x_1)} \n "
        print(info)
        file_name = "Lip" + self.base_name + ".txt"
        with open(file_name, "a") as file:
            file.write(info)
        return E

    def __init__(self, lower_bounds, upper_bounds, g, base_name, error=0.0001, N=1):

        self.base_name = base_name
        self.g = g
        self.bnds = [(lower_bounds[i], upper_bounds[i]) for i in range(len(lower_bounds))]
        self.bnds = self.bnds + self.bnds
        self.error = error
        self.N = N * len(self.bnds) * 1000  # default for Powell's mimimize

        # error=0.000001
        # N=10


# finding_Lipschitz(lower_bounds, upper_bounds, g, base_name)

# print(slope(g, [-2.34279779, - 2.51941202, - 2.56941876],
# [-2.44407475, - 2.27628564, - 2.55924496]))

# print(g([-2.34279779, - 2.51941202, - 2.56941876]))

# print(g([-2.44407475, - 2.27628564, - 2.55924496]))


# [2.55040391 2.36068774 0.28470773 2.75396691 2.36613446 0.29615437]
# -36.94207095524288


# def F(rect):
#     return F_K(rect, g, K)
# run_CMGDB(subdiv_min, subdiv_max, lower_bounds, upper_bounds, phase_periodic, F, base_name)
