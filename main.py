import fem.material as material
import fem.solver as solver
from fem.io.reader import S2KReader
from fem.io.writer import TextWriter
import math
import time
import numpy as np


mat = material.Steel('mat1', 2e8, 0.3, 8)
d = solver.Domain()
d.add_material(mat)

reader = S2KReader()

start_time = time.time()
reader.read('./data/input/sapinput.s2k', d)
print("---- %s reading seconds ----" % (time.time() - start_time))

start_time = time.time()
d.build()
print("---- %s building seconds ----" % (time.time() - start_time))

d.solve_linear()


writer = TextWriter()
writer.write('./data/output/deneme.txt', d)


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