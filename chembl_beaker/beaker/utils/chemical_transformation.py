__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem.AllChem import Compute2DCoords
from rdkit.Chem.AllChem import GenerateDepictionMatching2DStructure
from rdkit.Chem import Kekulize
from rdkit.Chem import SanitizeMol
from rdkit.Chem import AddHs
from rdkit.Chem import RemoveHs
from rdkit.Chem import GetSSSR
from rdkit.Chem import GetSymmSSSR
from rdkit.Chem.rdmolops import SanitizeFlags as sf
SANITIZE_ALL = sf.SANITIZE_ALL

# ----------------------------------------------------------------------------------------------------------------------


def _computeCoords(mol, force=False):
    if force or (not mol.GetNumConformers() or mol.GetConformer().Is3D()):
        Compute2DCoords(mol)
    return mol

# ----------------------------------------------------------------------------------------------------------------------


def _atomMapNumber(mol):
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return mol

# ----------------------------------------------------------------------------------------------------------------------


def _kekulize(mol):
    return Kekulize(mol)

# ----------------------------------------------------------------------------------------------------------------------


def _sanitize(mol, sanitizeOps=SANITIZE_ALL):
    return SanitizeMol(mol, sanitizeOps=sanitizeOps)

# ----------------------------------------------------------------------------------------------------------------------


def _addHs(mol, explicitOnly=False, addCoords=False):
    return AddHs(mol, explicitOnly=explicitOnly, addCoords=addCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _removeHs(mol, implicitOnly=False):
    return RemoveHs(mol, implicitOnly=implicitOnly)

# ----------------------------------------------------------------------------------------------------------------------


def _sssr(mol):
    return GetSSSR(mol)

# ----------------------------------------------------------------------------------------------------------------------


def _symmsssr(mol):
    return GetSymmSSSR(mol)

# ----------------------------------------------------------------------------------------------------------------------


def _getSubstructMatch(mol, patt):
    if mol.HasSubstructMatch(patt):
        return mol.GetSubstructMatch(patt)
    return []

# ----------------------------------------------------------------------------------------------------------------------


def _align(mols, pattern):
    GenerateDepictionMatching2DStructure(mols, pattern)
    return mols

# ----------------------------------------------------------------------------------------------------------------------

