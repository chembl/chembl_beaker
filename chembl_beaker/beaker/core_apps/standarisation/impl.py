__author__ = 'mnowotka'

from rdkit import Chem
from standardiser.break_bonds import apply as break_bonds
from standardiser.neutralise import apply as neutralise
from standardiser.rules import apply as rules
from standardiser.unsalt import is_nonorganic, is_salt
from standardiser.standardise import apply as standardise
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

def _break_bonds(data):
    mols = _parseMolData(data)
    res = _apply(mols, break_bonds)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _neutralise(data, balance):
    mols = _parseMolData(data)
    res = _apply(mols, neutralise, balance)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _rules(data):
    mols = _parseMolData(data)
    res = _apply(mols, rules)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _unsalt(data):
    mols = _parseMolData(data)
    res = _apply(mols, unsalt)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------

def _standardise(data):
    mols = _parseMolData(data)
    res = _apply(mols, standardise)
    return _getSDFString(res)

#-----------------------------------------------------------------------------------------------------------------------