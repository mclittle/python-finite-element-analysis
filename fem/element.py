import numpy as np
from numpy.linalg import *


class Solid8(object):

    def __init__(self, name, material, *nodes):
        self.name = name
        self.material = material
        self.nodes = nodes
        self.code = []
        self.NODOF = sum(n.NODOF for n in self.nodes)
        self.map = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]

    def __str__(self):
        s = 'Element{0} | '.format(self.name)
        return s

    @staticmethod
    def fi(r, s, t):
        return np.asarray([0.125 * (1 + item[0] * r) * (1 + item[1] * s) * (1 + item[2] * t) for item in map], np.float64)

    @staticmethod
    def df(r, s, t):
        map = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
        dFdr = np.asarray([0.125 * (item[0]) * (1 + item[1] * s) * (1 + item[2] * t) for item in map], np.float64)
        dFds = np.asarray([0.125 * (1 + item[0] * r) * (item[1]) * (1 + item[2] * t) for item in map], np.float64)
        dFdt = np.asarray([0.125 * (1 + item[0] * r) * (1 + item[1] * s) * (item[2]) for item in map], np.float64)
        return np.asarray([dFdr, dFds, dFdt]).T

    def get_jako(self, r, s, t):
        xv = np.asarray([self.nodes[i].coord for i in range(8)], np.float64).T
        return xv @ Solid8.df(r, s, t)

    def dFx(self, r, s, t):
        jako = self.get_jako(r, s, t)
        return Solid8.df(r, s, t) @ inv(jako)

    def get_c(self):
        E = self.material.e
        p = self.material.p
        G = self.material.g
        return inv(np.asarray([[1.0 / E, -p / E, -p / E, 0, 0, 0],
                               [-p / E, 1.0 / E, -p / E, 0, 0, 0],
                               [-p / E, -p / E, 1 / E, 0, 0, 0],
                               [0, 0, 0, 1.0 / G, 0, 0],
                               [0, 0, 0, 0, 1.0 / G, 0],
                               [0, 0, 0, 0, 0, 1.0 / G]], np.float64))

    def get_b(self, r, s, t):
        dfx = self.dFx(r, s, t)
        temp = np.zeros((6, 24))
        temp[0, 0:8] = dfx[:, 0]
        temp[1, 8:16] = dfx[:, 1]
        temp[2, 16:] = dfx[:, 2]

        temp[3, 0:8] = dfx[:, 1]
        temp[3, 8:16] = dfx[:, 0]

        temp[4, 0:8] = dfx[:, 2]
        temp[4, 16:] = dfx[:, 0]

        temp[5, 8:16] = dfx[:, 2]
        temp[5, 16:] = dfx[:, 1]

        return temp

    def get_dk(self, r, s, t):
        b = self.get_b(r, s, t)
        return b.T @ self.get_c() @ b

    def get_k(self):
        #map = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
        qp = np.asarray(self.map, np.float64) / 3.0**0.5
        temp = np.zeros((24, 24))
        for rst in qp:
            temp += self.get_dk(*rst) * det(self.get_jako(*rst))
        return temp

    def get_m(self):
        m = np.zeros((len(self.nodes), len(self.nodes)))
        #map = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]

        for rst in self.map:
            for i in range(len(self.nodes)):
                m[i, i] += self.material.ro * det(self.get_jako(*rst)) / len(self.nodes)

        return m

    # Get Displacement of nodes
    def get_d(self):
        d = np.zeros(self.NODOF)

        c = 0
        for n in self.nodes:
            for i, e in enumerate(n.q):
                d[i * len(self.nodes) + c] = n.u[i]
            c += 1
        return d

    # Get Sigma of Element
    def get_sigma(self, r, s, t):
        return self.get_c() @ self.get_b(r, s, t) @ self.get_d()

    # Apply element stresses to nodes
    def add_sigma_to_nodes(self):
        sigmas = [0] * len(self.nodes)
        for i, node in enumerate(self.nodes):
            sigmas[i] = self.get_sigma(self.map[i][0] / 3 ** 0.5, self.map[i][1] / 3 ** 0.5, self.map[i][2] / 3 ** 0.5)

        for i in range(0, 4):
            next = i + 4

            map = self.map

            x1 = map[i][0] / 3 ** 0.5
            y1 = map[i][1] / 3 ** 0.5
            z1 = map[i][2] / 3 ** 0.5
            x2 = map[next][0] / 3 ** 0.5
            y2 = map[next][1] / 3 ** 0.5
            z2 = map[next][2] / 3 ** 0.5

            rx1 = 0
            ry1 = 0
            rz1 = 0
            rx2 = 0
            ry2 = 0
            rz2 = 0

            if abs(x2 - x1) > 0:
                rx1 = (map[i][0] - x1) / (x2 - x1)
                rx2 = (map[next][0] - x1) / (x2 - x1)
            if abs(y2 - y1) > 0:
                ry1 = (map[i][1] - y1) / (y2 - y1)
                ry2 = (map[next][1] - y1) / (y2 - y1)
            if abs(z2 - z1) > 0:
                rz1 = (map[i][2] - z1) / (z2 - z1)
                rz2 = (map[next][2] - z1) / (z2 - z1)

            sigma1 = (rx1 + ry1 + rz1) * (sigmas[next] - sigmas[i]) + sigmas[i]
            sigma2 = (rx2 + ry2 + rz2) * (sigmas[next] - sigmas[i]) + sigmas[i]

            self.nodes[i].sigma += sigma1
            self.nodes[next].sigma += sigma2

'''class Solid8Nonlinear(Solid8):
    def get_jakoInit(self, r, s, t):
        xv = np.asarray([self.nodes[i].InitialCoord for i in range(8)], np.float64).T
        return xv @ Solid8.df(r, s, t)

    def get_F(self, r, s, t):
        return self.get_jako(r, s, t) @ inv(self.get_jakoInit(r, s, t))

    def get_sigma(self, r, s, t):
        lamb = (self.material.p * self.material.e) / ((1 + self.material.p) * (1 - 2 * self.material.p))
        nu = self.material.e / (2 * (1 + self.material.p))
        F = self.get_F(r, s, t)
        detF = det(F)
        I = np.ones([6, 6])
        return lamb @ (np.log(detF) / detF * I) + nu * (F @ F.T - I) / detF

    def get_dk(self, r, s, t):
        b = self.get_b(r, s, t)
        sigma = self.get_sigma(r, s, t)'''