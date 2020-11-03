import re
from fem.node import *
from fem.element import *


class S2KReader:

    def read(self, file_path, domain):
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
                    self.__read_joint_coordinates(line, domain)
                elif reading == 'CONNECTIVITY - SOLID':
                    self.__read_connectivity_solid(line, domain)
                elif reading == 'JOINT RESTRAINT ASSIGNMENTS':
                    self.__read_joint_restraint_assignments(line, domain)
                elif reading == 'JOINT LOADS - FORCE':
                    self.__read_joint_loads_force(line, domain)

    def __read_joint_coordinates(self, line, domain):
        name = re.search('Joint=(\w+)', line).group(1)

        x = float(re.search('GlobalX=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
        y = float(re.search('GlobalY=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
        z = float(re.search('GlobalZ=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))

        domain.nodes[name] = Node3(x, y, z)
        domain.nodes[name].name = name

    def __read_connectivity_solid(self, line, domain):
        name = re.search('Solid=(\w+)', line).group(1)
        node_list = []
        for match in re.findall('Joint[^=]+=(\w+)', line):
            node_list.append(domain.nodes[match])

        node_list = [node_list[0], node_list[1], node_list[3], node_list[2], node_list[4], node_list[5],
                     node_list[7], node_list[6]]

        domain.elements[name] = Solid8(name, domain.materials[0], *node_list)

    def __read_joint_restraint_assignments(self, line, domain):
        name = re.search('Joint=(\w+)', line).group(1)
        u1 = 1 if re.search('U1=(\w+)', line).group(1) == 'Yes' else 0
        u2 = 1 if re.search('U2=(\w+)', line).group(1) == 'Yes' else 0
        u3 = 1 if re.search('U3=(\w+)', line).group(1) == 'Yes' else 0

        domain.nodes[name].r = [u1, u2, u3]

    def __read_joint_loads_force(self, line, domain):
        name = re.search('Joint=(\w+)', line).group(1)
        f1 = float(re.search('F1=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
        f2 = float(re.search('F2=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
        f3 = float(re.search('F3=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))

        domain.nodes[name].p = [f1, f2, f3]
