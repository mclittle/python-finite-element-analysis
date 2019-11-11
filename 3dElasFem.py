import numpy as np
import re

from numpy.linalg import inv, det
np.set_printoptions(precision=5)

#Dictionary Keywords
code,rest, nodes, dof, Force, coord= 'code','rest', 'nodes', 'dof', 'Force', 'coord'

#Configuration Parameters
NODOF = 3  #Number of degree of freedom per node
NON   = 8  #Number of nodes
EDIM  = NODOF * NON
numbering = 'bydof' #'bynode'


MAT = {0: {'E': 2e8, 'p': 0.3}}

#Node and Element Dictionaries

ND = {}
ED = {}

# Read input file and create ND and ED dictionaries for algorithm
file = open('sapinput.s2k', 'r')
reading = ''
for line in file:
    if reading != '':
        if line.strip('\r\n') == '':
            reading = ''
            continue

        if reading == 'JOINT COORDINATES':
            name = re.search('Joint=(\w+)', line).group(1)
            x = float(re.search('GlobalX=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
            y = float(re.search('GlobalY=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
            z = float(re.search('GlobalZ=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
            ND[name] = {coord: [x, y, z]}

        elif reading == 'CONNECTIVITY - SOLID':
            name = re.search('Solid=(\w+)', line).group(1)
            node_list = []
            for match in re.findall('Joint[^=]+=(\w+)', line):
                node_list.append(match)

            node_list = [node_list[0], node_list[1], node_list[3], node_list[2], node_list[4], node_list[5], node_list[7], node_list[6]]

            ED[name] = {nodes: node_list, 'mat': MAT[0]}

        elif reading == 'JOINT RESTRAINT ASSIGNMENTS':
            name = re.search('Joint=(\w+)', line).group(1)
            u1 = 1 if re.search('U1=(\w+)', line).group(1) == 'Yes' else 0
            u2 = 1 if re.search('U2=(\w+)', line).group(1) == 'Yes' else 0
            u3 = 1 if re.search('U3=(\w+)', line).group(1) == 'Yes' else 0

            ND[name][rest] = [u1, u2, u3]

        elif reading == 'JOINT LOADS - FORCE':
            name = re.search('Joint=(\w+)', line).group(1)
            f1 = float(re.search('F1=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
            f2 = float(re.search('F2=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))
            f3 = float(re.search('F3=([-+]?[0-9]*[.,]?[0-9]+([eE][-+]?[0-9]+)?)', line).group(1).replace(',', '.'))

            ND[name][Force] = [f1, f2, f3]
        else:
            reading = ''
    else:
        parts = line.split(':')
        if parts[0] == 'TABLE':
            reading = parts[1].strip(' "\r\n')
            continue

for n in ND.values():
    print(n)
for e in ED.values():
    print(e)
'''
ND = {'0': {coord: [0, 0, 0], rest: [1, 1, 1]},
      '1': {coord: [2.5, 0, 0], rest: [1, 1, 1]},
      '2': {coord: [2.5, 3, 0], rest: [1, 1, 1]},
      '3': {coord: [0, 3, 0], rest: [1, 1, 1]},
      '4': {coord: [0, 0, 4], Force:[10, 0,0]},
      '5': {coord: [2.5, 0, 4]},
      '6': {coord: [2.5, 3, 4]},
      '7': {coord: [0, 3, 4]}
      }

for i in range(1, 10):
    ND[str(i * 4 + 4)] = {coord: [0, 0, i * 4.0 + 4.0], Force: [1000, 0, 0]}
    ND[str(i * 4 + 5)] = {coord: [2.5, 0, i * 4.0 + 4.0]}
    ND[str(i * 4 + 6)] = {coord: [2.5, 3, i * 4.0 + 4.0]}
    ND[str(i * 4 + 7)] = {coord: [0, 3, i * 4.0 + 4.0]}

ED = {'0': {'mat': MAT[0], nodes: ['0', '1', '2', '3', '4', '5', '6', '7']}
      }
for i in range(1, 10):
    ED[str(i)] = {'mat':MAT[0], nodes:[str(i*4), str(i*4 + 1), str(i*4 + 2), str(i*4 + 3),
                                       str(i*4 + 4), str(i*4 + 5), str(i*4 + 6), str(i*4 + 7)]}
'''
#*********** 3D Elasticity Element *********************


def XV(e):return np.asarray([ND[e[nodes][i]][coord] for i in range(8)], float).T

map = [[-1,-1,-1],[1, -1,-1],[1, 1, -1], [-1, 1,-1], [-1, -1, 1], [1, -1, 1], [1, 1, 1],[-1, 1, 1]]


def F(r, s, t): return np.asarray([0.125*(1+item[0]*r)*(1+item[1]*s)*(1+item[2]*t) for item in map], float)
def dFdr(r, s, t): return np.asarray([0.125*(item[0])*(1+item[1]*s)*(1+item[2]*t) for item in map], float)
def dFds(r, s, t): return np.asarray([0.125*(1+item[0]*r)*(item[1])*(1+item[2]*t) for item in map], float)
def dFdt(r, s, t): return np.asarray([0.125*(1+item[0]*r)*(1+item[1]*s)*(item[2]) for item in map], float)
def dF(r, s, t): return np.asarray([dFdr(r,s,t), dFds(r,s,t), dFdt(r,s,t)]).T
def JAKO(e, rst): return XV(e) @ dF(*rst)
def INVJAKO(e, rst): return inv(JAKO(e, rst))
def dFx(e, rst): return dF(*rst) @ INVJAKO(e, rst)

def C(e):
    E = e['mat']['E']
    p = e['mat']['p']
    G = 0.5 * E / (1 + p)
    return inv(np.asarray([[1.0/E, -p/E,   -p/E, 0, 0 ,0],
                           [-p/E , 1.0/E,  -p/E, 0, 0, 0],
                           [-p/E,  -p/E,   1/E,  0, 0, 0],
                           [0, 0, 0, 1.0/G, 0, 0 ],
                           [0, 0, 0, 0, 1.0 / G, 0],
                           [0, 0, 0, 0, 0, 1.0 / G]], float))

def B(e, rst):
    dfx = dFx(e, rst)
    temp = np.zeros((6, 24))
    temp[0,0:8]=dfx[:, 0]
    temp[1,8:16] = dfx[:, 1]
    temp[2,16:] = dfx[:, 2]

    temp[3, 0:8] = dfx[:, 1]
    temp[3, 8:16] = dfx[:, 0]

    temp[4, 0:8] = dfx[:, 2]
    temp[4, 16:] = dfx[:, 0]

    temp[5, 8:16] = dfx[:, 2]
    temp[5, 16:] = dfx[:, 1]

    return temp

def dK(e, rs):
    b = B(e, rs)
    return b.T @ C(e) @ b

qp = np.asarray(map, float) / 3.0**0.5

def KG(e):
    temp = np.zeros((24,24))
    for rst in qp:
        temp += dK(e, rst) * det(JAKO(e, rst))
    return temp

def sigma_xyz(e, rst):
    return C(e) @ B(e, rst) @ LHS[e[code]]


#****************  PREPARE  **************************************
#Apply missing values to node dictionary (ND)
for n in ND.values():
    n[dof] = [0]*NODOF
    n[Force] = n.get(Force, [0]*NODOF)
    n[rest] = n.get(rest, [0]*NODOF)
#Enumerating Nodal degree-of-freedoms
M = 0
N = 0
for kere in [0, 1]:
    for n in ND.values():
        for i, r in enumerate(n[rest]):
            if r==kere: n[dof][i]=M; M+=1
    if kere==0: N = M
#Finding code vectors

print(ND)

for e in ED.values():
    e[code] = []
    if numbering=='bydof':
        for i in range(NODOF):
            for n in e[nodes]: e[code] += [ND[n][dof][i]]
    else:
        for n in e[nodes]: e[code] += ND[n][dof]


#Print some info
print('Number of displacement unknowns: ', N)
print('Number of total unknowns: ', M)


#********************** ASSEMBLE ********************
#Create ZERO vectors and matrices
RHS = np.zeros(M)
LHS = np.zeros(M)
A = np.zeros(M**2).reshape(M,M)
#*************************************************************
#Burasi A yi oluşruruyor Onumuzdeki hafta anlatilacak        *
#*************************************************************
for e in ED.values():
    columns = np.asarray(e[code]*EDIM)
    rows = columns.reshape(EDIM, EDIM).T.reshape(EDIM**2)
    A[rows, columns] += KG(e).reshape(EDIM**2)
#*************************************************************
#Constructing the RHS
for n in ND.values():
    for i, d in enumerate(n[dof]): RHS[d] = n[Force][i]

#****************** SOLVE  ********************
#Solving the Linear System equation
LHS[0:N] = np.linalg.solve(A[0:N,0:N], RHS[0:N])
RHS[N:M] = A[N:M,:]@LHS

#***************** OUTPUT *****************************
#print('yerdeğiştirmeler: ',LHS)
#print('mesnet tepkileri: ',RHS)

print('========= Yer Değiştirmeler =========')
for key, value in ND.items():
    s = 'node{' + key + '} | '
    for q in value['dof']:
        s += str(LHS[q]) + ' | '
    print(s)

#print(sigma_xyz(ED['23'], [1, 1, 1]))


# Create OUTPUT FILE #
file = open('output.txt', 'w')
output = ''

output += 'START:NODES \n'
for key, value in ND.items():
    output += key + ';'
    for c in value[coord]:
        output += str(c) + ';'
    output += '\n'
output += 'END \n'

output += 'START:ELEMENTS \n'
for key, value in ED.items():
    output += key + ';'
    for n in value[nodes]:
        output += n + ';'
    output += '\n'
output += 'END \n'

output += 'START:DISPLACEMENTS \n'
for key, value in ND.items():
    output += key + ';'
    for q in value['dof']:
        output += str(LHS[q]) + ';'
    output += '\n'
output += 'END \n'

file.write(output)

#map2 = [[-1,-1,-1],[1, -1,-1],[1, 1, -1], [-1, 1,-1], [-1, -1, 1], [1, -1, 1], [1, 1, 1],[-1, 1, 1]]
#def F2(r, s, t): return np.asarray([0.125*(1+item[0]*r)*(1+item[1]*s)*(1+item[2]*t) for item in map2], float)

#print(F2(0, 0, 0))