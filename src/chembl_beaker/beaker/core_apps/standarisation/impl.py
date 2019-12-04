__author__ = 'efelix'

from rdkit import Chem
from chembl_structure_pipeline import standardizer, checker
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData


#-----------------------------------------------------------------------------------------------------------------------

def _check(data, loadMol=False, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, checker.check_molblock)
    return str(res)

#-----------------------------------------------------------------------------------------------------------------------

def _get_parent(data, loadMol=False, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, standardizer.get_parent_molblock)
    return str(res)

#-----------------------------------------------------------------------------------------------------------------------

def _standardise(data, loadMol=False, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, standardizer.standardize_molblock)
    return res

#-----------------------------------------------------------------------------------------------------------------------
