#/usr/bin/env pypy
class Bond:
    def __init__(self, bond):
        self.atoms = [bond[0], bond[1]]
        self.bond_type = bond[2]

    @property
    def bonded_atoms(self):
        return self.atoms

    @property
    def bond_type(self):
        return self.bond_type
