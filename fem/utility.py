import math

def xfrange(start, stop, step):
    i = 0.0
    factor = 1 / step
    while start + i * step <= stop:
        print(factor)
        yield math.ceil((start + i * step)*factor)/factor
        i += 1.0