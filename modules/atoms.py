class Atom:
    def __init__(self, atom):
        self._symbol = atom[0]
        self._index = atom[1]
        self._x = atom[2]
        self._y = atom[3]
        self._z = atom[4]

    @property
    def symbol(self):
        return self._symbol

    @property
    def index(self):
        return self._index

    @property
    def position(self):
        return self._x, self._y, self._z
