__author__ = 'mnowotka'

from rdkit import Chem
from standardiser.break_bonds import run as break_bonds
from standardiser.neutralise import run as neutralise
from standardiser.rules import run as rules
from standardiser.unsalt import is_nonorganic, is_salt
from standardiser.standardise import run as standardise
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString

#-----------------------------------------------------------------------------------------------------------------------

def unsalt(mol):
    non_salt_frags = []
    for n, frag in enumerate(Chem.GetMolFrags(mol, asMols=True), 1):
        if is_nonorganic(frag) or is_salt(frag):
            continue
        non_salt_frags.append(frag)

    if len(non_salt_frags) != 1:
        return []

    return non_salt_frags[0]

#-----------------------------------------------------------------------------------------------------------------------

def _break_bonds(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = _apply(mols, break_bonds)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _neutralise(data, balance, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = _apply(mols, neutralise, balance)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _rules(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = _apply(mols, rules)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _unsalt(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = _apply(mols, unsalt)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _standardise(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = _apply(mols, standardise)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------
