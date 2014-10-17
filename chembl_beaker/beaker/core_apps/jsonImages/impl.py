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

def _ctab2json(data, size, legend):
    return _mols2json(_parseMolData(data), size,legend)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2json(data, size, legend):
    return _mols2json(_parseSMILESData(data), size,legend)

#-----------------------------------------------------------------------------------------------------------------------
