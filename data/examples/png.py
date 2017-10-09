from rdkit import Chem
from rdkit.Chem import AllChem


suppl = Chem.SDMolSupplier("2.sdf")
ms = [x for x in suppl if x is not None]
for m in ms: tmp=AllChem.Compute2DCoords(m)
from rdkit.Chem import Draw
for x in range(len(ms)):
	Draw.MolToFile(ms[x],str(x)+'.png',subImgSize=(500,500))


