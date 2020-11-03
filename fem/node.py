
import numpy as np


class Node3(object):
    NODOF = 3

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.InitialCoord = [x, y, z]
        self.coord = [x, y, z]
        self.r = [0]*3
        self.p = [0]*3
        self.q = [0]*3
        self.u = [0]*3
        self.sigma = [0] * 6

    def __str__(self):
        return 'Node{0} ({1}, {2}, {3}) | {4} | {5}'.format(self.name, self.x, self.y, self.z, self.r, self.p)





