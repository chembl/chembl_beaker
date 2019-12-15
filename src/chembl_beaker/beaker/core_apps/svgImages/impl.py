__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

try:
    import cairo
    cffi = False
except ImportError:
    import cairocffi
    cairocffi.install_as_pycairo()
    cffi = True
    import io
    import cairo
    if not hasattr(cairo, 'HAS_PDF_SURFACE'):
        cairo.HAS_PDF_SURFACE = False
    if not hasattr(cairo, 'HAS_SVG_SURFACE'):
        cairo.HAS_SVG_SURFACE = True


from itertools import cycle, islice
from rdkit import Chem
from beaker.utils.functional import _apply, _call
from beaker.utils.io import _parseMolData, _parseSMILESData
from beaker.utils.chemical_transformation import _computeCoords, _atomMapNumber, _kekulize
from beaker.utils.io import _getMatches
from rdkit.Chem.Draw import rdMolDraw2D
import beaker.utils.chemical_transformation as ct


# ----------------------------------------------------------------------------------------------------------------------


def _mols2svg(mols, size, kekulize=False, atomMapNumber=False, computeCoords=False, highlightAtomLists=None):

    if not mols:
        return ''
    else:
        # process only one molecule
        mols = mols[0:1]

    _call(mols, 'UpdatePropertyCache', strict=False)
    _apply(mols, ct._sssr)

    if computeCoords:
        _apply(mols, _computeCoords, True)
    if atomMapNumber:
        _apply(mols, _atomMapNumber)

    if kekulize:
        _apply(mols, _kekulize)

    if highlightAtomLists:
        highlightAtomLists = highlightAtomLists[0]

    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    drawer.DrawMolecule(mols[0], highlightAtoms=highlightAtomLists)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2svg(data, size, loadMol=True, useRDKitChemistry=False, kekulize=False,
              atomMapNumber=False, computeCoords=False):
    return _mols2svg(_parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry),
        size, kekulize=kekulize, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2svg(data, size, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                titleLine=True, sanitize=True, kekulize=True, atomMapNumber=False):
    return _mols2svg(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize), size,
        kekulize=kekulize, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _inchi2svg(inchis, size, kekulize=True, atomMapNumber=False,
               computeCoords=False):
    mols = _apply(inchis.split(), Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _mols2svg(mols, size, kekulize=kekulize,
                        atomMapNumber=atomMapNumber, computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightSmilesFragmentSVG(data, smarts, size, computeCoords=False, delimiter=' ', smilesColumn=0,
                                nameColumn=1, titleLine=True, sanitize=True, kekulize=True,
                                atomMapNumber=False, force=False):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                            smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    matches = _getMatches(mols, smarts, force)
    return _mols2svg(mols, size, kekulize=kekulize,
                     atomMapNumber=atomMapNumber, computeCoords=computeCoords, highlightAtomLists=matches)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightCtabFragmentSVG(data, smarts, size, loadMol=True, useRDKitChemistry=False,
                              kekulize=False, atomMapNumber=False, computeCoords=False,
                              force=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    matches = _getMatches(mols, smarts, force)
    return _mols2svg(mols, size, kekulize=kekulize,
                     atomMapNumber=atomMapNumber, computeCoords=computeCoords, highlightAtomLists=matches)

# ----------------------------------------------------------------------------------------------------------------------
