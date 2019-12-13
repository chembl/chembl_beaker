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
from rdkit.Chem.Draw import rdMolDraw2D, MolDrawOptions
import beaker.utils.chemical_transformation as ct


# ----------------------------------------------------------------------------------------------------------------------


def _mols2svg(mols, size, legend, kekulize=False, fitImage=True, atomMapNumber=False,
              computeCoords=False, highlightAtomLists=None):

    if not mols:
        return ''

    _call(mols, 'UpdatePropertyCache', strict=False)
    _apply(mols, ct._sssr)

    if computeCoords:
        _apply(mols, _computeCoords, True)
    if atomMapNumber:
        _apply(mols, _atomMapNumber)

    labels = [x for x in islice(cycle(legend), len(mols))] if isinstance(legend, (list, tuple)) else \
             [x for x in islice(cycle([legend]), len(mols))]
    legends = [x.GetProp("_Name") if (x.HasProp("_Name") and x.GetProp("_Name")) else labels[idx]
               for idx, x in enumerate(mols)]

    molsPerRow = min(len(mols), 4)
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow:
        nRows += 1

    if kekulize:
        _apply(mols, _kekulize)

    highlightBondLists = []
    if highlightAtomLists:
        for mol, highlightAtomList in zip(mols, highlightAtomLists):
            highlightBondList = []
            for bnd in mol.GetBonds():
                if bnd.GetBeginAtomIdx() in highlightAtomList and bnd.GetEndAtomIdx() in highlightAtomList:
                    highlightBondList.append(bnd.GetIdx())
            highlightBondLists.append(highlightBondList)
    if not highlightBondLists:
        highlightBondLists = None

    panelx = size
    panely = size
    canvasx = panelx * molsPerRow
    canvasy = panely * nRows
    drawer = rdMolDraw2D.MolDraw2DSVG(canvasx, canvasy, panelx, panely)
    drawer.DrawMolecules(mols, highlightAtoms=highlightAtomLists, highlightBonds=highlightBondLists, legends=legends)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2svg(data, size, legend, loadMol=True, useRDKitChemistry=False, kekulize=False,
              fitImage=True, atomMapNumber=False, computeCoords=False):
    return _mols2svg(_parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry),
        size, legend, kekulize=kekulize, fitImage=fitImage, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2svg(data, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                titleLine=True, sanitize=True, kekulize=True, fitImage=True, atomMapNumber=False):
    return _mols2svg(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize), size, legend,
        kekulize=kekulize, fitImage=fitImage, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _inchi2svg(inchis, size, legend, kekulize=True, fitImage=True, atomMapNumber=False,
               computeCoords=False):
    mols = _apply(inchis.split(), Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _mols2svg(mols, size, legend, kekulize=kekulize, fitImage=fitImage,
                        atomMapNumber=atomMapNumber, computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightSmilesFragmentSVG(data, smarts, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0,
                                nameColumn=1, titleLine=True, sanitize=True, kekulize=True,
                                fitImage=True, atomMapNumber=False, force=False):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                            smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    matches = _getMatches(mols, smarts, force)
    return _mols2svg(mols, size, legend, kekulize=kekulize, fitImage=fitImage,
                     atomMapNumber=atomMapNumber, computeCoords=computeCoords, highlightAtomLists=matches)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightCtabFragmentSVG(data, smarts, size, legend, loadMol=True, useRDKitChemistry=False,
                              kekulize=False, fitImage=True, atomMapNumber=False, computeCoords=False,
                              force=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    matches = _getMatches(mols, smarts, force)
    return _mols2svg(mols, size, legend, kekulize=kekulize, fitImage=fitImage,
                     atomMapNumber=atomMapNumber, computeCoords=computeCoords, highlightAtomLists=matches)

# ----------------------------------------------------------------------------------------------------------------------
