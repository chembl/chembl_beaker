__author__ = 'efelix'

from rdkit import Chem
from chembl_structure_pipeline import standardizer, checker
from beaker.utils.functional import _apply
from beaker.utils.io import _parseMolData
import json


#-----------------------------------------------------------------------------------------------------------------------

def _check(data, loadMol=False, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, checker.check_molblock)
    return json.dumps(res)

#-----------------------------------------------------------------------------------------------------------------------

def _get_parent(data, loadMol=False, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, standardizer.get_parent_molblock)
    res_dict = {k: v for e in res for k, v in zip(('parent_molblock', 'exclude'), e)}
    return json.dumps(res_dict)

#-----------------------------------------------------------------------------------------------------------------------

def _standardise(data, loadMol=False, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, standardizer.standardize_molblock)
    return res

#-----------------------------------------------------------------------------------------------------------------------
