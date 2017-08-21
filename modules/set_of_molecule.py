import rdkit
from rdkit import Chem
from molecules import Molecule
from sys import exit


class Set_of_molecule:
    def __init__(self, file):
        with open(file, 'r') as if_sdf:
            if if_sdf.readlines()[-1].split(None, 1)[0] != "$$$$":
                exit("ERROR! File is not sdf file!")
        set_of_molecules_rdkit = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        list_with_molecules = []
        for molecule_rdkit in set_of_molecules_rdkit:
            rdkit.Chem.Kekulize(molecule_rdkit)
            mol_data = {}
            list_with_atoms = []
            number_of_atoms = molecule_rdkit.GetNumAtoms()
            for i in range(0, number_of_atoms):
                pos = molecule_rdkit.GetConformer().GetAtomPosition(i)
                list_with_atoms.append((molecule_rdkit.GetAtomWithIdx(i).GetSymbol(), i+1, pos.x, pos.y, pos.z))
            mol_data["atoms"] = list_with_atoms
            list_with_bonds = []
            for i in range(0, number_of_atoms):
                for j in range(i, number_of_atoms):
                    try:
                        list_with_bonds.append((i+1, j+1, str(molecule_rdkit.GetBondBetweenAtoms(i, j).GetBondType())))
                    except:
                        pass
            mol_data["bonds"] = list_with_bonds
            mol_data["total_charge"] = rdkit.Chem.GetFormalCharge(molecule_rdkit)
            list_with_molecules.append(Molecule((molecule_rdkit.GetProp('_Name'), mol_data)))
        self.list_with_molecules = list_with_molecules
        self.number_of_molecules = len(list_with_molecules)

    def molecules(self, returned_molecule):
        if returned_molecule > self.number_of_molecules:
            exit("ERROR! You want to more molecules from set, then molecules are in set!")
        return self.list_with_molecules[:returned_molecule]

    @property
    def number_of_molecules(self):
        return self.number_of_molecules


"""
format of molecule in self.list_with_molecules

('NSC_102026',
 {'atoms': [('C', 1, -0.6392, -0.5604, 0.0011),
            ('C', 2, -0.0547, 1.6192, -0.001),
            ('C', 3, 0.4564, 0.2999, -0.0),
            ('O', 4, 2.1798, -1.2657, 0.0021),
            ('C', 5, 1.8718, -0.0884, 0.0005),
            ('N', 6, -1.7472, 0.2083, 0.0009),
            ('N', 7, -1.3558, 1.5513, -0.0004),
            ('N', 8, 2.8305, 0.8589, -0.0008),
            ('N', 9, -0.5984, -1.9399, -0.0028),
            ('C', 10, -3.1281, -0.2808, 0.0011),
            ('H', 11, 0.5322, 2.5258, -0.0014),
            ('H', 12, 2.5851, 1.7973, -0.002),
            ('H', 13, 3.7659, 0.6022, -0.0004),
            ('H', 14, -1.4238, -2.4495, -0.0058),
            ('H', 15, 0.2556, -2.3999, -0.0026),
            ('H', 16, -3.4706, -0.4017, 1.0288),
            ('H', 17, -3.7677, 0.4367, -0.513),
            ('H', 18, -3.1735, -1.2412, -0.5124)],
  'bonds': [(1, 3, 'DOUBLE'),
            (1, 6, 'SINGLE'),
            (1, 9, 'SINGLE'),
            (2, 3, 'SINGLE'),
            (2, 7, 'DOUBLE'),
            (2, 11, 'SINGLE'),
            (3, 5, 'SINGLE'),
            (4, 5, 'DOUBLE'),
            (5, 8, 'SINGLE'),
            (6, 7, 'SINGLE'),
            (6, 10, 'SINGLE'),
            (8, 12, 'SINGLE'),
            (8, 13, 'SINGLE'),
            (9, 14, 'SINGLE'),
            (9, 15, 'SINGLE'),
            (10, 16, 'SINGLE'),
            (10, 17, 'SINGLE'),
            (10, 18, 'SINGLE')],
  'total_charge': 0})

"""
