__author__ = 'mnowotka'

from rdkit import Chem
from rdkit.Chem.rdmolops import SanitizeFlags as sf
SANITIZE_ALL = sf.SANITIZE_ALL
from chembl_beaker.beaker.utils.functional import _apply, _call
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString

#-----------------------------------------------------------------------------------------------------------------------

def _atomRings(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _call(_call(mols, 'GetRingInfo'), 'AtomRings')

#-----------------------------------------------------------------------------------------------------------------------

def _bondRings(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _call(_call(mols, 'GetRingInfo'), 'BondRings')

#-----------------------------------------------------------------------------------------------------------------------

def _isAtomInRing(data, index, size, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _call(_call(mols, 'GetRingInfo'), 'IsAtomInRingOfSize', index, size)

#-----------------------------------------------------------------------------------------------------------------------

def _isBondInRing(data, index, size, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _call(_call(mols, 'GetRingInfo'), 'IsBondInRingOfSize', index, size)

#-----------------------------------------------------------------------------------------------------------------------

def _numAtomRings(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    ring_infos = _call(mols, 'GetRingInfo')
    return [[ring_info.NumAtomRings(atom.GetIdx()) for atom in mol.GetAtoms()] for (mol, ring_info) in zip(mols, ring_infos)]

#-----------------------------------------------------------------------------------------------------------------------

def _numBondRings(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    ring_infos = _call(mols, 'GetRingInfo')
    return [[ring_info.NumBondRings(bond.GetIdx()) for bond in mol.GetBonds()] for (mol, ring_info) in zip(mols, ring_infos)]

#-----------------------------------------------------------------------------------------------------------------------

def _numRings(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _call(_call(mols, 'GetRingInfo'), 'NumRings')

#-----------------------------------------------------------------------------------------------------------------------
