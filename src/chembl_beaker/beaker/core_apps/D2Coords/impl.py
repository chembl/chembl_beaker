__author__ = 'efelix'

from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString, _parseSMILESData


# ----------------------------------------------------------------------------------------------------------------------


def _2D22D(mol):
    Chem.rdDepictor.SetPreferCoordGen(True)
    Chem.rdDepictor.Compute2DCoords(mol)
    return mol

# ----------------------------------------------------------------------------------------------------------------------


def _ctab22D(data, multi, sanitize=False, removeHs=False, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    optimisedMols = _apply(mols, _2D22D)
    return _getSDFString(optimisedMols)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles22D(data, multi, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
               titleLine=True, sanitize=False):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter, smilesColumn=smilesColumn,
        nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    optimisedMols = _apply(mols, _2D22D)
    return _getSDFString(optimisedMols)

# ----------------------------------------------------------------------------------------------------------------------
