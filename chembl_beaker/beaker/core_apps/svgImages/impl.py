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
from chembl_beaker.beaker.draw import cairoCanvas
from chembl_beaker.beaker import draw
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords

#-----------------------------------------------------------------------------------------------------------------------

def _mols2svg(mols,size,legend):

    _apply(mols, _computeCoords)

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
        leg = mol.GetProp("_Name") if mol.HasProp("_Name") else legend
        draw.MolToImage(mol, size=(size,size), legend=leg, canvas=canv)
        canv.flush()
        ctx.translate(-tx, -ty)
    surf.finish()
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _ctab2svg(data,size,legend):
    return _mols2svg(_parseMolData(data), size, legend)

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2svg(data,size,legend):
    return _mols2svg(_parseSMILESData(data), size, legend)

#-----------------------------------------------------------------------------------------------------------------------
