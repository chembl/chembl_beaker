__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import StringIO
from rdkit.Chem import SDMolSupplier
from rdkit.Chem import SmilesMolSupplier
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import SDWriter
from rdkit.Chem import SmilesWriter
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords

#-----------------------------------------------------------------------------------------------------------------------

def _parseMolData(data):
    suppl = SDMolSupplier()
    suppl.SetData(str(data))
    return [x for x in suppl if x]

#-----------------------------------------------------------------------------------------------------------------------

def _parseSMILESData(data, computeCoords=False):
    suppl = SmilesMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl if x]
    if not mols:
        mols = [MolFromSmiles(data)]
    if computeCoords:
        _apply(mols, _computeCoords, True)
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

def _getSMILESStream(f, mols):
    w = SmilesWriter(f, isomericSmiles=True)
    for mol in mols:
        w.write(mol)
    w.flush()

#-----------------------------------------------------------------------------------------------------------------------

def _getSMILESString(mols):
    sio = StringIO.StringIO()
    _getSMILESStream(sio, mols)
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------
