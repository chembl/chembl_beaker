__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import StringIO
from itertools import cycle, islice
from chembl_beaker.beaker import draw
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords, _atomMapNumber
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData

#-----------------------------------------------------------------------------------------------------------------------

def _mols2imageStream(mols, f, format, size, legend):
    labels = [x for x in islice(cycle(legend), len(mols))] if isinstance(legend, (list, tuple)) else \
             [x for x in islice(cycle([legend]), len(mols))]
    legends = [x.GetProp("_Name") if (x.HasProp("_Name") and x.GetProp("_Name")) else labels[idx]
               for idx, x in enumerate(mols)]
    image = draw.MolsToGridImage(mols, molsPerRow=min(len(mols), 4), subImgSize=(size,size), legends=legends)
    image.save(f, format)

#-----------------------------------------------------------------------------------------------------------------------

def _mols2imageString(mols, size, legend, format, atomMapNumber=False, computeCoords=False):
    if not mols:
        return ''
    if computeCoords:
        _apply(mols, _computeCoords, True)
    if atomMapNumber:
        _apply(mols, _atomMapNumber)
    imageData = StringIO.StringIO()
    _mols2imageStream(mols, imageData, format, size, legend)
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _ctab2image(data, size, legend, sanitize=True, removeHs=True, strictParsing=True, atomMapNumber=False,
                computeCoords=False):
    return _mols2imageString(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
        size, legend, 'PNG', atomMapNumber, computeCoords)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2image(data, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                  titleLine=True, sanitize=True, atomMapNumber=False):
    return _mols2imageString(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize), size, legend, 'PNG',
        atomMapNumber, computeCoords)

#-----------------------------------------------------------------------------------------------------------------------
