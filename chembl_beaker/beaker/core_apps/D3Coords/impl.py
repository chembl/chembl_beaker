__author__ = 'mnowotka'

from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString, _parseSMILESData

#-----------------------------------------------------------------------------------------------------------------------

def _2D23D(mol, multi, mmff=False):
    molH = Chem.AddHs(mol)
    if mmff:
        optiFunc = AllChem.MMFFOptimizeMolecule
    else:
        optiFunc = AllChem.UFFOptimizeMolecule
    if multi:
        confIds=AllChem.EmbedMultipleConfs(molH, multi)
        for confId in confIds:
          optiFunc(molH, confId=confId)
    else:
        AllChem.EmbedMolecule(molH)
        optiFunc(molH)

    return molH

#-----------------------------------------------------------------------------------------------------------------------

def _ctab23D(data, multi, mmff, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    optimisedMols = _apply(mols, _2D23D, multi, mmff)
    return _getSDFString(optimisedMols)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles23D(data, multi, mmff, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
               titleLine=True, sanitize=True):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter, smilesColumn=smilesColumn,
        nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    optimisedMols = _apply(mols, _2D23D, multi, mmff)
    return _getSDFString(optimisedMols)

#-----------------------------------------------------------------------------------------------------------------------
