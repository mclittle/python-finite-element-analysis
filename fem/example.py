
def __create_domain():
    d = solver.Domain()

    mat = material.Steel('mat1', 2e8, 0.3, 8)
    d.add_material(mat)

    return d


def one_element():
    domain = __create_domain()

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
        domain.nodes[n.name] = n

    elm1 = Solid8('1', domain.materials[0], n1, n2, n3, n4, n5, n6, n7, n8)

    domain.elements['elm1'] = elm1

    return domain


def two_element():
    domain = __create_domain()

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
        domain.nodes[n.name] = n

    elm1 = Solid8('1', self.materials[0], n1, n2, n3, n4, n5, n6, n7, n8)
    elm2 = Solid8('2', self.materials[0], n5, n6, n7, n8, n9, n10, n11, n12)

    domain.elements['elm1'] = elm1
    domain.elements['elm2'] = elm2

    return domain
