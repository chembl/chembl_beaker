__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem.AllChem import Compute2DCoords

#-----------------------------------------------------------------------------------------------------------------------

def _computeCoords(mol, force=False):
    if force or (not mol.GetNumConformers() or mol.GetConformer().Is3D()):
        Compute2DCoords(mol)
    return mol

#-----------------------------------------------------------------------------------------------------------------------

def _atomMapNumber(mol):
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return mol

#-----------------------------------------------------------------------------------------------------------------------
