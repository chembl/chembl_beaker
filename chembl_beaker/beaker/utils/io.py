__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
import StringIO
from rdkit.Chem import SDMolSupplier
from rdkit.Chem import SmilesMolSupplier
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import SDWriter
from rdkit.Chem import SmilesWriter
from rdkit.Chem import SanitizeMol
from rdkit.Chem import SanitizeFlags as sf
from chembl_beaker.beaker.utils.functional import _apply, _call
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords

#-----------------------------------------------------------------------------------------------------------------------

def _parseFlag(data):
    try:
        return bool(int(data))
    except ValueError:
        if isinstance(data, basestring):
            l = data.lower()
            if l == 't' or l == 'true':
                return True
            if l == 'f' or l == 'false':
                return False
        if not data:
            return False
        return None

#-----------------------------------------------------------------------------------------------------------------------

def _parseMolData(data, sanitize=True, removeHs=True, strictParsing=True):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = SDMolSupplier(fpath, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = [x for x in suppl if x]
    os.remove(fpath)
    return res

#-----------------------------------------------------------------------------------------------------------------------

def _parseSMILESData(data, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                     sanitize=True):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = SmilesMolSupplier(fpath, delimiter=delimiter, smilesColumn=smilesColumn, nameColumn=nameColumn,
        titleLine=titleLine, sanitize=sanitize)
    mols = [x for x in suppl if x]
#    if not mols:
#        mols = [MolFromSmiles(data, sanitize=sanitize)]
    if computeCoords:
        _apply(mols, _computeCoords, True)
    _call(mols, 'SetProp', '_Name', '')
    os.remove(fpath)
    return mols

#-----------------------------------------------------------------------------------------------------------------------

def _getSDFStream(f, mols):
    w = SDWriter(f)
    for m in mols:
        w.write(m)
    w.flush()

#-----------------------------------------------------------------------------------------------------------------------

def _getSDFString(mols):
    sio = StringIO.StringIO()
    _getSDFStream(sio, mols)
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _getSMILESStream(f, mols, delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                     kekuleSmiles=False):
    w = SmilesWriter(f, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
        isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)
    for mol in mols:
        w.write(mol)
    w.flush()

#-----------------------------------------------------------------------------------------------------------------------

def _getSMILESString(mols, delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                     kekuleSmiles=False):
    sio = StringIO.StringIO()
    _getSMILESStream(sio, mols, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
        isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------
