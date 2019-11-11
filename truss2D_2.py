import numpy as np

ND = {0:{'x':[0,0], 'q':[0,11], 'r':[0,1], 'force': [0, 0]},
      1:{'x':[2000,0], 'q':[1,2], 'r':[0,0], 'force':[0, -15]},
      2:{'x':[4000,0], 'q':[3,12], 'r':[0,1], 'force': [0, 0]},
      3:{'x':[6000,0], 'q':[4,5], 'r':[0,0], 'force':[0,-10]},
      4:{'x':[1000,1500], 'q':[13,6], 'r':[1,0], 'force': [0, 0]},
      5:{'x':[3000,1500], 'q':[7,8], 'r':[0,0], 'force': [0, 0]},
      6:{'x':[5000,1500], 'q':[9,10], 'r':[0,0], 'force':[5, 0]}
      }

# E = 200000 A = 0.25x0.50
ED = {0:{'conn':[0,1], 'EA':25000000000},
      1:{'conn':[1,2], 'EA':25000000000},
      2:{'conn':[2,3], 'EA':25000000000},
      3:{'conn':[0,4], 'EA':25000000000},
      4:{'conn':[4,1], 'EA':25000000000},
      5:{'conn':[1,5], 'EA':25000000000},
      6:{'conn':[5,2], 'EA':25000000000},
      7:{'conn':[2,6], 'EA':25000000000},
      8:{'conn':[6,3], 'EA':25000000000},
      9:{'conn':[4,5], 'EA':25000000000},
      10:{'conn':[5,6], 'EA':25000000000}
      }

'Wrapper' 'Utility'
def x(n):
    return ND[n]['x']

def conn(e):
    return ED[e]['conn']

def xij(e, ij):
    return x(conn(e)[ij])
    #return ND[ED[e]['conn'][ij]]['x']

def xi(e):
    return xij(e, 0)

def xj(e):
    return xij(e, 1)


def kod(e):
      left_node_id = ED[e]['conn'][0]
      rigth_node_id = ED[e]['conn'][1]
      return ND[left_node_id]['q'] + ND[rigth_node_id]['q']


#The following functions are your homework to implement
#I have already implemented the last one for you :))

def L(e):
    return ((xj(e)[0]-xi(e)[0])**2+(xj(e)[1]-xi(e)[1])**2)**0.5

def nx(e):
      return (xj(e)[0]-xi(e)[0])/L(e)

def ny(e):
      return (xj(e)[1]-xi(e)[1])/L(e)

#print('nx:', nx(1))
#print('ny:', ny(1))

def T(e):
    '''Transformation matrix of the element with id=e'''
    a = nx(e)
    b = ny(e)
    return np.asarray([ [a,  b, 0, 0],
             [-b, a, 0, 0],
             [0, 0,  a, b],
             [0, 0, -b, a] ])

def KL(e):
    EA_L = ED[e]['EA']/L(e)
    return np.asarray([[EA_L, 0, -EA_L, 0],
            [0, 0, 0, 0],
            [-EA_L, 0, EA_L, 0],
            [0, 0, 0, 0]])


def KG(e):
    return T(e).T @ KL(e) @ T(e)

#[]  List
#{} Dict
#() Tuple

A = np.zeros((14, 14))

for e in range(10):
    a = KG(e).flatten()
    cv = kod(e)
    cols = cv*4
    rows = np.asarray(cv*4).reshape(4,4).T.flatten().tolist()
    A[rows, cols] += a


RHS = np.zeros(14)
for n in ND.values():
    for i, d in enumerate(n['q']): RHS[d] = n['force'][i]

u = np.zeros(14)
u[0:11] =np.linalg.inv(A[0:11,0:11])@RHS[0:11]
print(RHS)
R = A@u
print(u)
print(R)










