#/usr/bin/env pypy
class Atom:
    def __init__(self, atom):
        self.symbol = atom[0]
        self.index = atom[1]
        self.x = atom[2]
        self.y = atom[3]
        self.z = atom[4]

    @property
    def symbol(self):
        return self.symbol

    @property
    def index(self):
        return self.index

    @property
    def position(self):
        return self.x, self.y, self.z
