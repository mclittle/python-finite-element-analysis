import re
import numpy as np
from numpy.linalg import *
from FEM.Node import *
from FEM.Element import *
from FEM.Utility import *
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

    def from_sap(self, file_path):
        file = open(file_path, 'r')
        reading = ''
        for line in file:
            parts = line.split(':')
            if parts[0] == 'TABLE':
                reading = parts[1].strip(' "\r\n')
                continue

            if reading != '':
                if not line.strip():
                    reading = ''
                    continue

                if reading == 'JOINT COORDINATES':
                        name = re.search('Joint=(\w+)', line).group(1)
                        x = float(
                            re.search('GlobalX=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',',
                                                                                                                   '.'))
                        y = float(
                            re.search('GlobalY=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',',
                                                                                                                   '.'))
                        z = float(
                            re.search('GlobalZ=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',',
                                                                                                                   '.'))
                        self.nodes[name] = Node3(x, y, z)
                        self.nodes[name].name = name

                elif reading == 'CONNECTIVITY - SOLID':
                    name = re.search('Solid=(\w+)', line).group(1)
                    node_list = []
                    for match in re.findall('Joint[^=]+=(\w+)', line):
                        node_list.append(self.nodes[match])

                    node_list = [node_list[0], node_list[1], node_list[3], node_list[2], node_list[4], node_list[5],
                                 node_list[7], node_list[6]]
                    self.elements[name] = Solid8(name, self.materials[0], *node_list)

                    #print(self.elements[name], node_list[0], node_list[1], node_list[2], node_list[3], node_list[4], node_list[5], node_list[6], node_list[7])

                elif reading == 'JOINT RESTRAINT ASSIGNMENTS':
                    name = re.search('Joint=(\w+)', line).group(1)
                    u1 = 1 if re.search('U1=(\w+)', line).group(1) == 'Yes' else 0
                    u2 = 1 if re.search('U2=(\w+)', line).group(1) == 'Yes' else 0
                    u3 = 1 if re.search('U3=(\w+)', line).group(1) == 'Yes' else 0

                    self.nodes[name].r = [u1, u2, u3]

                elif reading == 'JOINT LOADS - FORCE':
                    name = re.search('Joint=(\w+)', line).group(1)
                    f1 = float(
                        re.search('F1=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
                    f2 = float(
                        re.search('F2=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
                    f3 = float(
                        re.search('F3=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))

                    self.nodes[name].p = [f1, f2, f3]

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

    def get_kg(self):
        edim = 24
        kg = np.zeros(self.m ** 2).reshape(self.m, self.m)
        for e in self.elements.values():
            columns = np.asarray(e.code * edim)
            rows = columns.reshape(edim, edim).T.reshape(edim ** 2)
            kg[rows, columns] += e.get_k().reshape(edim ** 2)

        return kg

    def get_p(self):
        p = np.zeros(self.m)
        for n in self.nodes.values():
            for i, d in enumerate(n.q): p[d] = n.p[i]

        return p

    def get_m(self):
        m = np.zeros((self.m, self.m))
        for e in self.elements.values():
            elm_m = e.get_m()
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
        kg = self.get_kg()
        p = self.get_p()
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

    #def solve_nonlinear(self):


    def example1(self):
        n1 = Node3(0, 0, 0)
        n1.name = 'n1'
        n1.r = [1, 1, 1]
        n2 = Node3(10, 0, 0)
        n2.name = 'n2'
        n2.r = [1, 1, 1]
        n3 = Node3(10, 10, 0)
        n3.name = 'n3'
        n3.r = [1, 1, 1]
        n4 = Node3(0, 10, 0)
        n4.name = 'n4'
        n4.r = [1, 1, 1]
        n5 = Node3(0, 0, 10)
        n5.name = 'n5'
        n6 = Node3(10, 0, 10)
        n6.name = 'n6'
        n7 = Node3(10, 10, 10)
        n7.name = 'n7'
        n8 = Node3(0, 10, 10)
        n8.name = 'n8'
        n9 = Node3(0, 0, 20)
        n9.name = 'n9'
        n9.p = [0, 0, 1]
        n10 = Node3(10, 0, 20)
        n10.name = 'n1'
        n10.p = [0, 0, 1]
        n11 = Node3(10, 10, 20)
        n11.name = 'n11'
        n11.p = [0, 0, 1]
        n12 = Node3(0, 10, 20)
        n12.name = 'n12'
        n12.p = [0, 0, 1]
        nodes = [n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12]
        for n in nodes:
            self.nodes[n.name] = n

        elm1 = Solid8('1', self.materials[0], n1, n2, n3, n4, n5, n6, n7, n8)
        elm2 = Solid8('2', self.materials[0], n5, n6, n7, n8, n9, n10, n11, n12)

        self.elements['elm1'] = elm1
        self.elements['elm2'] = elm2

    def one_element(self):
        n1 = Node3(0, 0, 0)
        n1.r = [1, 1, 1]
        n2 = Node3(0.25, 0, 0)
        n2.r = [1, 1, 1]
        n3 = Node3(0.25, 0.25, 0)
        n3.r = [1, 1, 1]
        n4 = Node3(0, 0.25, 0)
        n4.r = [1, 1, 1]
        n5 = Node3(0, 0, 0.25)
        n5.p = [0, 0, 5]
        n6 = Node3(0.25, 0, 0.25)
        n6.p = [0, 0, 5]
        n7 = Node3(0.25, 0.25, 0.25)
        n7.p = [0, 0, 5]
        n8 = Node3(0, 0.25, 0.25)
        n8.p = [0, 0, 5]
        n1.name = '1'
        n2.name = '2'
        n3.name = '3'
        n4.name = '4'
        n5.name = '5'
        n6.name = '6'
        n7.name = '7'
        n8.name = '8'

        nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
        for n in nodes:
            self.nodes[n.name] = n

        elm1 = Solid8('1', self.materials[0], n1, n2, n3, n4, n5, n6, n7, n8)

        self.elements['elm1'] = elm1

    def solve_explicit(self, dt, t0, t1, p):
        k = self.get_kg()
        m = self.get_m()
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
        k = self.get_kg()
        m = self.get_m()
        uprev = np.zeros(self.m)
        u = np.zeros(self.m)
        v = np.zeros(self.m)
        a = np.zeros(self.m)
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
        k = self.get_kg()
        m = self.get_m()
        u = np.zeros(self.m)
        v = np.zeros(self.m)
        a = np.zeros(self.m)
        output = self.write_to_string()
        max_times = 100
        inv_m = inv(m[0:self.n, 0:self.n])
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
        m = self.get_m()
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

    def write_to_file(self, outputfile):
        # Create OUTPUT FILE #
        file = open(outputfile, 'w')
        output = ''

        output += 'START:NODES \n'
        for key, n in self.nodes.items():
            output += key + ';'
            for c in n.coord:
                output += str(c) + ';'
            for c in n.r:
                output += str(c) + ';'
            for c in n.p:
                output += str(c) + ';'
            output += '\n'
        output += 'END:NODES \n'

        output += 'START:ELEMENTS \n'
        for key, elm in self.elements.items():
            output += key + ';'
            for n in elm.nodes:
                output += n.name + ';'
            output += '\n'
        output += 'END:ELEMENTS \n'

        output += 'START:DISPLACEMENTS \n'
        for n in self.nodes.values():
            output += n.name + ';'
            for u in n.u:
                output += str(u) + ';'
            output += '\n'
        output += 'END:DISPLACEMENTS \n'


        output += 'START:STRESSES \n'
        for n in self.nodes.values():
            output += n.name + ';'
            for u in n.sigma:
                output += str(u) + ';'
            output += '\n'
        output += 'END:STRESSES \n'

        file.write(output)

    def write_to_file_forJS(self, outputfile):
        file = open(outputfile, 'w')
        output = ''

        output += 'START:VERTICES \n'
        for key, n in self.nodes.items():
            output += key + ';'
            for e, c in enumerate(n.coord):
                output += str(c + n.u[e]) + ';'
            output += str(n.u[2])
            output += '\n'
        output += 'END:VERTICES \n'

        output += 'START:ELEMENTS \n'
        for key, elm in self.elements.items():
            output += key + ';'
            for n in elm.nodes:
                output += n.name + ';'
            output += '\n'
        output += 'END:ELEMENTS \n'

        file.write(output)