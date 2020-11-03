

class TextWriter:

    def __init__(self):
        self.stream = None
        self.domain = None

    def write(self, filepath, domain):
        self.stream = open(filepath, 'w')
        self.domain = domain

        self.__write_nodes()
        self.__writer_elements()
        self.__writer_displacements()
        self.__writer_stresses()

        self.stream.close()

    def __write_nodes(self):
        self.stream.write('START:NODES \n')
        for key, n in self.domain.nodes.items():
            line = key + ';'
            for c in n.coord:
                line += str(c) + ';'
            for c in n.r:
                line += str(c) + ';'
            for c in n.p:
                line += str(c) + ';'
            line += '\n'
            self.stream.write(line)
        self.stream.write('END:NODES \n')

    def __writer_elements(self):
        self.stream.write('START:ELEMENTS \n')
        for key, elm in self.domain.elements.items():
            line = key + ';'
            for n in elm.nodes:
                line += n.name + ';'
            line += '\n'
            self.stream.write(line)
        self.stream.write('END:ELEMENTS \n')

    def __writer_displacements(self):
        self.stream.write('START:DISPLACEMENTS \n')
        for n in self.domain.nodes.values():
            line = n.name + ';'
            for u in n.u:
                line += str(u) + ';'
            line += '\n'
            self.stream.write(line)
        self.stream.write('END:DISPLACEMENTS \n')

    def __writer_stresses(self):
        self.stream.write('START:STRESSES \n')
        for n in self.domain.nodes.values():
            line = n.name + ';'
            for u in n.sigma:
                line += str(u) + ';'
            line += '\n'
            self.stream.write(line)
        self.stream.write('END:STRESSES \n')