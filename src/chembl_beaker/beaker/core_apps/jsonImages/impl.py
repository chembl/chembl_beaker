__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker.core_apps.jsonImages.jsonCanvas import MolToJSON
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords

#-----------------------------------------------------------------------------------------------------------------------

def _mol2json(mol, size, legend):
    leg = mol.GetProp("_Name") if mol.HasProp("_Name") else legend
    return MolToJSON(_computeCoords(mol), size=(size,size), legend=leg)

#-----------------------------------------------------------------------------------------------------------------------

def _mols2json(mols,size,legend):
    return _apply(mols, _mol2json, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

def _ctab2json(data, size, legend, sanitize=True, removeHs=True, strictParsing=True):
    return _mols2json(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
        size,legend)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2json(data, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                 titleLine=True, sanitize=True):
    return _mols2json(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize), size,legend)

#-----------------------------------------------------------------------------------------------------------------------
