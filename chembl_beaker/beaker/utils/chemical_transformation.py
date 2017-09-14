__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker.utils.functional import _apply, _call
from rdkit.Chem.AllChem import Compute2DCoords
from rdkit.Chem.AllChem import GenerateDepictionMatching2DStructure
from rdkit.Chem import Kekulize
from rdkit.Chem import AddHs
from rdkit.Chem import RemoveHs
from rdkit.Chem import GetSSSR
from rdkit.Chem import GetSymmSSSR
from rdkit.Chem import SanitizeMol
from rdkit.Chem import AdjustQueryProperties
from rdkit.Chem import AdjustQueryParameters
from rdkit.Chem.rdmolops import SanitizeFlags as sf
from itertools import compress
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


def _adjustQuery(pattern):
    params = AdjustQueryParameters()
    params.adjustDegree = False
    params.makeBondsGeneric = True
    return AdjustQueryProperties(pattern, params)

# ----------------------------------------------------------------------------------------------------------------------


def _getSubstructMatch(mol, patt, force=False):
    if mol.HasSubstructMatch(patt):
        return mol.GetSubstructMatch(patt)
    if not force:
        return []
    patt = _adjustQuery(patt)
    if mol.HasSubstructMatch(patt):
        return mol.GetSubstructMatch(patt)
    return []

# ----------------------------------------------------------------------------------------------------------------------


def _align(mols, pattern, force=False):

    if all(mol.HasSubstructMatch(pattern) for mol in mols):
        _apply(mols, GenerateDepictionMatching2DStructure, pattern)
        return mols

    if not force:
        return

    pattern = _adjustQuery(pattern)

    matches = _call(mols, 'HasSubstructMatch', pattern)

    if all(matches):
        _apply(mols, GenerateDepictionMatching2DStructure, pattern)
        return mols

    if any(matches):
        _apply(mols, GenerateDepictionMatching2DStructure, pattern, acceptFailure=True)
        return compress(mols, matches)

# ----------------------------------------------------------------------------------------------------------------------

