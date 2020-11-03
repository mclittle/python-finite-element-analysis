
import numpy as np
from numpy.linalg import *
from fem.node import *
from fem.element import *
from fem.utility import *
import time
np.set_printoptions(16)


class Domain(object):
    nodes = {}
    elements = {}
    materials = []
    n = 0
    m = 0

    def add_material(self, material):
        self.materials.append(material)

    def build(self):
        for restrained in [0, 1]:
            for n in self.nodes.values():
                for i, r in enumerate(n.r):
                    if r == restrained:
                        n.q[i] = self.m
                        self.m += 1
            if restrained == 0:
                self.n = self.m

        for e in self.elements.values():
            for i in range(3):
                for n in e.nodes:
                    e.code += [n.q[i]]

    def __get_kg(self):
        edim = 24
        kg = np.zeros(self.m ** 2).reshape(self.m, self.m)
        for e in self.elements.values():
            columns = np.asarray(e.code * edim)
            rows = columns.reshape(edim, edim).T.reshape(edim ** 2)
            kg[rows, columns] += e.get_k().reshape(edim ** 2)

        return kg

    def __get_p(self):
        p = np.zeros(self.m)
        for n in self.nodes.values():
            for i, d in enumerate(n.q): p[d] = n.p[i]

        return p

    def __get_m(self):
        m = np.zeros((self.m, self.m))
        for e in self.elements.values():
            elm_m = e.__get_m()
            for n in e.nodes:
                for i, d in enumerate(n.q):
                    m[d, d] += elm_m[i, i]
        return m

    def calculate_sigmas(self):
        for elm in self.elements.values():
            elm.add_sigma_to_nodes()

        for node in self.nodes.values():
            elm_count = sum(1 for i in self.elements.values() if node in i.nodes)
            node.sigma = node.sigma / elm_count

    def solve_linear(self):
        result = np.zeros(self.m)

        start_time = time.time()
        kg = self.__get_kg()
        p = self.__get_p()
        print("---- %s assembly seconds ----" % (time.time() - start_time))

        start_time = time.time()
        result[0:self.n] = np.linalg.solve(kg[0:self.n, 0:self.n], p[0:self.n])
        print("---- %s solving seconds ----" % (time.time() - start_time))

        for n in self.nodes.values():
            for i, e in enumerate(n.q):
                n.u[i] = result[e]

    def __str__(self):
        s = 'Nodes \n'
        for n in self.nodes:
            s += n + '\n'
        s += 'Elements \n'
        for e in self.elements:
            s += e + '\n'
        return s

    def solve_explicit(self, dt, t0, t1, p):
        k = self.__get_kg()
        m = self.__get_m()
        u = np.zeros(self.m)
        v = np.zeros(self.m)
        output = self.write_to_string()
        max_times = (t1 - t0) / 100 / dt
        inv_m = inv(m[0:self.n, 0:self.n])

        file = open('output_explicit.txt', 'w')
        times = 0
        for t in xfrange(t0, t1, dt):
            pv = np.zeros(self.m)
            for n in self.nodes.values():
                for i, d in enumerate(n.q):
                    pv[d] = p(t) * n.p[i]

            #print(t)
            times += 1

            if times > max_times or t - dt < dt * 10.0 or t1 - t < dt:
                output += 'START:DISPLACEMENTS:t=' + str(t) + ' \n'
                for n in self.nodes.values():
                    output += n.name + ';'
                    output += ';'.join(str(q) for q in u[n.q])
                    output += '\n'
                output += 'END:DISPLACEMENTS \n'
                file.write(output)
                output = ''
                print(t)
                times = 0

            dv = (inv_m @ (pv[0:self.n] - k[0:self.n, 0:self.n] @ u[0:self.n])) * dt
            du = v[0:self.n] * dt
            u[0:self.n] += du
            v[0:self.n] += dv

        for n in self.nodes.values():
            n.u = u[n.q]

    def solve_explicit_central(self, dt, t0, t1, f):
        k = self.__get_kg()
        m = self.__get_m()
        u = np.zeros(self.m)
        output = self.write_to_string()
        max_times = (t1 - t0) / 100 / dt - 1

        p0 = np.zeros(self.m)
        for n in self.nodes.values():
            for i, d in enumerate(n.q):
                p0[d] = f(t0) * n.p[i]

        kkcap = inv(m / (dt ** 2))
        A = k - 2 * m / (dt ** 2)
        B = m / (dt ** 2)

        a0 = inv(m) @ (p0 - k @ u)

        uprev = u + a0 * dt * dt / 2

        file = open('output_central.txt', 'w')
        times = 0
        for t in xfrange(t0, t1, dt):
            fv = np.zeros(self.m)
            for n in self.nodes.values():
                for i, d in enumerate(n.q):
                    fv[d] = f(t)*n.p[i]

            fcap = fv[0:self.n] - A[0:self.n, 0:self.n] @ u[0:self.n] - B[0:self.n, 0:self.n] @ uprev[0:self.n]
            unext = kkcap[0:self.n, 0:self.n] @ fcap[0:self.n]
            uprev[0:self.n] = u[0:self.n]
            u[0:self.n] = unext

            if times > max_times or t - dt < dt * 10.0 or t1 - t < dt:
                output += 'START:DISPLACEMENTS:t=' + str(t) + ' \n'
                for n in self.nodes.values():
                    output += n.name + ';'
                    output += ';'.join(str(q) for q in u[n.q])
                    output += '\n'
                output += 'END:DISPLACEMENTS \n'
                file.write(output)
                output = ''
                print(t)
                if times > max_times:
                    times = 0

            times += 1

        for n in self.nodes.values():
            n.u = u[n.q]

    def solve_implicit(self, dt, ttimes, p):
        k = self.__get_kg()
        m = self.__get_m()
        u = np.zeros(self.m)
        v = np.zeros(self.m)
        a = np.zeros(self.m)
        output = self.write_to_string()
        max_times = 100
        beta = 0.25
        gamma = 0.5
        a1 = (1 / (beta * dt * dt)) * m[0:self.n, 0:self.n]
        a2 = (1 / (beta * dt)) * m[0:self.n, 0:self.n]
        a3 = (1 / (2 * beta) - 1) * m[0:self.n, 0:self.n]
        kcap = (k[0:self.n, 0:self.n] + a1)

        inv_kcap = inv(kcap)
        p0 = np.zeros(self.m)
        for n in self.nodes.values():
            for i, d in enumerate(n.q):
                p0[d] = p(0) * n.p[i]

        file = open('output_soliddynamic.txt', 'w')
        times = 0
        for tt in range(0, ttimes + 1):
            t = round(tt * dt, 10)
            pv = np.zeros(self.m)
            pnext = np.zeros(self.m)
            for n in self.nodes.values():
                for i, d in enumerate(n.q):
                    pv[d] = p(t) * n.p[i]
                    pnext[d] = p(t + dt) * n.p[i]

            pcap = pnext[0:self.n] + a1 @ u[0:self.n] + a2 @ v[0:self.n] + a3 @ a[0:self.n]
            unext = inv_kcap @ pcap
            vnext = gamma / (beta * dt) * (unext - u[0:self.n]) + (1 - gamma / beta) * v[0:self.n] + dt * (
                1 - gamma / (2 * beta)) * a[0:self.n]
            anext = 1 / (beta * dt * dt) * (unext - u[0:self.n]) - (1 / (beta * dt)) * v[0:self.n] - (1 / (
                2 * beta) - 1) * a[0:self.n]

            if times == 0:
                output += 'START:DISPLACEMENTS:t=' + str(t) + ' \n'
                for n in self.nodes.values():
                    output += n.name + ';'
                    output += ';'.join(str(q) for q in u[n.q])
                    output += '\n'
                output += 'END:DISPLACEMENTS \n'
                file.write(output)
                output = ''
                print(t)

            times += 1

            if times >= max_times:
                times = 0

            if t == 1 and True:
                print('t: ', t)
                print("u")
                print('2 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['2'].q]), '\\\\')
                print('3 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['3'].q]), '\\\\')
                print('6 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['6'].q]), '\\\\')
                print('8 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['8'].q]), '\\\\')
                print('10 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['10'].q]), '\\\\')
                print('12 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['12'].q]), '\\\\')
                print('87 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['87'].q]), '\\\\')
                print('88 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['88'].q]), '\\\\')
                print('99 & ', ' & '.join(str(round(i, 20)) for i in u[self.nodes['99'].q]), '\\\\')
                print("v")
                print('2 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['2'].q]), '\\\\')
                print('3 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['3'].q]), '\\\\')
                print('6 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['6'].q]), '\\\\')
                print('8 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['8'].q]), '\\\\')
                print('10 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['10'].q]), '\\\\')
                print('12 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['12'].q]), '\\\\')
                print('87 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['87'].q]), '\\\\')
                print('88 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['88'].q]), '\\\\')
                print('99 & ', ' & '.join(str(round(i, 20)) for i in v[self.nodes['99'].q]), '\\\\')
                print("a")
                print('2 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['2'].q]), '\\\\')
                print('3 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['3'].q]), '\\\\')
                print('6 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['6'].q]), '\\\\')
                print('8 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['8'].q]), '\\\\')
                print('10 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['10'].q]), '\\\\')
                print('12 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['12'].q]), '\\\\')
                print('87 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['87'].q]), '\\\\')
                print('88 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['88'].q]), '\\\\')
                print('99 & ', ' & '.join(str(round(i, 20)) for i in a[self.nodes['99'].q]), '\\\\')
                # exit()

            u[0:self.n] = unext
            v[0:self.n] = vnext
            a[0:self.n] = anext

        for n in self.nodes.values():
            n.u = u[n.q]

    def write_to_string(self):
        output = ''

        output += 'START:NODES \n'
        m = self.__get_m()
        for key, n in self.nodes.items():
            output += key + ';'
            for c in n.coord:
                output += str(c) + ';'
            for c in n.r:
                output += str(c) + ';'
            for c in n.p:
                output += str(c) + ';'
            for c in n.q:
                output += str(m[c, c]) + ';'
            output += '\n'
        output += 'END:NODES \n'

        output += 'START:ELEMENTS \n'
        for key, elm in self.elements.items():
            output += key + ';'
            for n in elm.nodes:
                output += n.name + ';'
            output += '\n'
        output += 'END:ELEMENTS \n'

        return output
