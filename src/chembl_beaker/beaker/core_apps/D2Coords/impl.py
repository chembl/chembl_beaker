__author__ = 'efelix'

from rdkit import Chem
from rdkit.Chem import AllChem
from beaker.utils.functional import _apply
from beaker.utils.io import _parseMolData, _getSDFString, _parseSMILESData
import json


# ----------------------------------------------------------------------------------------------------------------------


def _2D22D(mol):
    Chem.rdDepictor.SetPreferCoordGen(True)
    Chem.rdDepictor.Compute2DCoords(mol)
    return mol

# ----------------------------------------------------------------------------------------------------------------------


def _check3Dcoords(mol):
    conf = mol.GetConformer()
    if conf.Is3D():
        return True

    for i in range(mol.GetNumAtoms()):
        if abs(conf.GetAtomPosition(i).z) >= 0.0001:
            return True
    return False

# ----------------------------------------------------------------------------------------------------------------------


def _ctab22D(data, loadMol=True, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=True, useRDKitChemistry=False)
    optimisedMols = _apply(mols, _2D22D)
    return _getSDFString(optimisedMols)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles22D(data, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
               titleLine=True, sanitize=False):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter, smilesColumn=smilesColumn,
        nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    optimisedMols = _apply(mols, _2D22D)
    return _getSDFString(optimisedMols)

# ----------------------------------------------------------------------------------------------------------------------


def _is3D(data, loadMol=True, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    flags = _apply(mols, _check3Dcoords)
    return json.dumps(flags)

# ----------------------------------------------------------------------------------------------------------------------
