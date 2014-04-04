__author__ = 'mnowotka'

from rdkit.Chem import AllChem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString, _parseSMILESData

#-----------------------------------------------------------------------------------------------------------------------

def _2D23D(mol, multi):
    if multi:
        confIds=AllChem.EmbedMultipleConfs(mol, multi)
        for confId in confIds:
          AllChem.UFFOptimizeMolecule(mol, confId=confId)
    else:
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

#-----------------------------------------------------------------------------------------------------------------------

def _ctab23D(data, multi):
    mols = _parseMolData(data)
    _apply(mols, _2D23D, multi)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles23D(data, multi):
    mols = _parseSMILESData(data)
    _apply(mols, _2D23D, multi)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------
