import FEM.Material as material
import FEM.Solver as solver
import math
import time
import numpy as np

mat = material.Steel('mat1', 2e8, 0.3, 8)

d = solver.Domain()
d.add_material(mat)
start_time = time.time()
d.from_sap('C:/Users/mclittle/Desktop/Programming/github/mclittle/sap_models/sapkiris2/sapinput.s2k')
print("---- %s reading seconds ----" % (time.time() - start_time))
#d.one_element()

start_time = time.time()
d.build()
print("---- %s building seconds ----" % (time.time() - start_time))

def p(t):
    t_max = 2
    t_min = 0
    count = 20
    interval = (t_max - t_min) / count
    t_start = math.floor(t / interval) * interval
    t_end = t_start + interval

    p1 = math.cos(t_start * math.pi)
    p2 = math.cos(t_end * math.pi)

    return (p2 - p1) * ((t-t_start) / (t_end - t_start)) + p1

#d.solve_implicit2()
#d.solve_implicit(1e-4, 10000, p)

d.solve_linear()
#d.calculate_sigmas()

#print(d.nodes['5'].sigma[0])

d.write_to_file('./data/outputs/output_solidstatic.txt')
# d.write_to_file_forJS('outputJS.txt')
