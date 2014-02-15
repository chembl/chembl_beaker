__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData, _getSMILESString, _getSDFString
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords

#-----------------------------------------------------------------------------------------------------------------------

def _ctab2smiles(data):
    return _getSMILESString(_parseMolData(data))

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2ctab(data):
    return _getSDFString(_parseSMILESData(data, True))

#-----------------------------------------------------------------------------------------------------------------------

def _inchi2ctab(inchis):
    mols = _apply(inchis.split(),Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------

def _ctab2inchi(data):
    return '\n'.join(_apply(_parseMolData(data), Chem.MolToInchi))

#-----------------------------------------------------------------------------------------------------------------------

def _inchi2inchiKey(inchis):
    return '\n'.join([Chem.InchiToInchiKey(inch) for inch in inchis.split()])

#-----------------------------------------------------------------------------------------------------------------------
