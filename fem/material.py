

class Steel(object):
    def __init__(self, name, elasticity_modulus, poisson_ratio, unit_mass):
        self.name = name
        self.e = elasticity_modulus
        self.p = poisson_ratio
        self.ro = unit_mass
        self.g = 0.5 * self.e / (1 + self.p)