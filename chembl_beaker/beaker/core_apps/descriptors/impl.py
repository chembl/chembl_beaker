__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import Descriptors
from chembl_beaker.beaker.utils.functional import _call, _apply
from chembl_beaker.beaker.utils.io import _parseMolData

#-----------------------------------------------------------------------------------------------------------------------

def _desc(mol, name):
    if name and isinstance(name, basestring) and hasattr(Descriptors, name):
        return getattr(Descriptors, name)(mol)

#-----------------------------------------------------------------------------------------------------------------------

def _desc_list(mol, names):
    descriptors = dict()
    for name,fn in Descriptors.descList:
        if not names or name in names:
            descriptors[name]=fn(mol)
    return descriptors

#-----------------------------------------------------------------------------------------------------------------------

def _getNumAtoms(data):
    return _call(_parseMolData(data), "GetNumAtoms")

#-----------------------------------------------------------------------------------------------------------------------

def _getLogP(data):
    return _apply(_parseMolData(data), _desc, 'MolLogP')

#-----------------------------------------------------------------------------------------------------------------------

def _getTPSA(data):
    return _apply(_parseMolData(data), _desc, 'TPSA')

#-----------------------------------------------------------------------------------------------------------------------

def _getMolWt(data):
    return _apply(_parseMolData(data), _desc, 'MolWt')

#-----------------------------------------------------------------------------------------------------------------------

def _getDescriptors(data, params={}):
    ds=params.get('descrs','')
    if ds!='':
        ds = ds.split(',')
    return _apply(_parseMolData(data), _desc_list, ds)

#-----------------------------------------------------------------------------------------------------------------------
