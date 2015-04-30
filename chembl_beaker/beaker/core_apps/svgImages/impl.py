__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

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


import StringIO
from rdkit import Chem
from chembl_beaker.beaker.draw import cairoCanvas
from chembl_beaker.beaker import draw
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords, _atomMapNumber

#-----------------------------------------------------------------------------------------------------------------------

def _mols2svg(mols, size, legend, kekulize=True, wedgeBonds=True, fitImage=True, atomMapNumber=False, computeCoords=False):

    if computeCoords:
        _apply(mols, _computeCoords, True)

    if atomMapNumber:
        _apply(mols, _atomMapNumber)

    molsPerRow=min(len(mols),4)
    totalWidth=molsPerRow*size
    totalHeight=molsPerRow*size
    if cffi and cairocffi.version <= (1,10,0) :
        imageData = io.BytesIO()
    else:
        imageData = StringIO.StringIO()
    surf = cairo.SVGSurface(imageData,totalWidth,totalHeight)
    ctx = cairo.Context(surf)
    for i, mol in enumerate(mols):
        tx = size*(i%molsPerRow)
        ty = size*(i//molsPerRow)
        ctx.translate(tx, ty)
        canv = cairoCanvas.Canvas(ctx=ctx, size=(size,size), imageType='svg')
        leg = mol.GetProp("_Name") if (mol.HasProp("_Name") and mol.GetProp("_Name")) else legend
        draw.MolToImage(mol, size=(size,size), legend=leg, canvas=canv, kekulize=kekulize, wedgeBonds=wedgeBonds,
               fitImage=fitImage)
        canv.flush()
        ctx.translate(-tx, -ty)
    surf.finish()
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _ctab2svg(data, size, legend, sanitize=True, removeHs=True, strictParsing=True, kekulize=True, wedgeBonds=True,
              fitImage=True, atomMapNumber=False, computeCoords=False):
    return _mols2svg(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
        size, legend, kekulize=kekulize, wedgeBonds=wedgeBonds, fitImage=fitImage, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2svg(data, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                titleLine=True, sanitize=True, kekulize=True, wedgeBonds=True, fitImage=True, atomMapNumber=False):
    return _mols2svg(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize), size, legend,
        kekulize=kekulize, wedgeBonds=wedgeBonds, fitImage=fitImage, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

#-----------------------------------------------------------------------------------------------------------------------

def _inchi2svg(inchis,size,legend, kekulize=True, wedgeBonds=True, fitImage=True, atomMapNumber=False, computeCoords=False):
    mols = _apply(inchis.split(), Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _mols2svg(mols, size, legend, kekulize=kekulize, wedgeBonds=wedgeBonds, fitImage=fitImage,
                        atomMapNumber=atomMapNumber, computeCoords=computeCoords)

#-----------------------------------------------------------------------------------------------------------------------