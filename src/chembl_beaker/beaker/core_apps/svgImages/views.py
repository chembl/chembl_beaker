__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from beaker import app
from bottle import request, response
from beaker.utils.io import _parseFlag
from beaker.core_apps.svgImages.impl import _ctab2svg, _smiles2svg, _inchi2svg
from beaker.core_apps.svgImages.impl import _highlightCtabFragmentSVG, _highlightSmilesFragmentSVG
from beaker.core_apps.svgImages.impl import _smiles2SimilarityMapSVG, _sdf2SimilarityMapSVG
from bottle import request, response

# ----------------------------------------------------------------------------------------------------------------------


def ctab2svgView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 300))
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', False))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))

    response.content_type = 'image/svg+xml'
    return _ctab2svg(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2svg', method=['OPTIONS', 'POST'], name="ctab2svg")
def ctab2svg():
    """
Converts CTAB to SVG vector graphic. CTAB is either single molfile or SDF file. Size is the optional size of
image in pixels (default value is 300 px). 
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "size=400" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@mcs.sdf" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs.sdf" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs_no_coords.sdf" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs.sdf" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs.sdf" -F "size=400" ${BEAKER_ROOT_URL}ctab2svg > out.svg    
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return ctab2svgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2svgView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 300))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/svg+xml'
    return _smiles2svg(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2svg', method=['OPTIONS', 'POST'], name="smiles2svg")
def smiles2svg():
    """
Converts SMILES to SVG vector graphic. This method accepts single or multiple SMILES or *.smi file. Size is the
optional size of image in pixels (default value is 300 px). 
cURL examples:

    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > out.svg    
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smiles2svgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def highlightSmilesFragmentSvgView(data, params):

    kwargs = dict()
    smarts = params.get('smarts', '')
    kwargs['size'] = int(params.get('size', 300))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['force'] = _parseFlag(params.get('force', True))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/svg+xml'
    return _highlightSmilesFragmentSVG(data, smarts, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/highlightSmilesFragmentSvg', method=['OPTIONS', 'POST'], name="highlightSmilesFragmentSvg")
def highlightSmilesFragmentSvg():
    """
Converts SMILES to SVG vector graphic with a highlighted fragment described as SMARTS. 
This method accepts SMARTS and single or multiple SMILES or *.smi file. Size is the
optional size of image in pixels (default value is 300 px). 
cURL examples:

    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@CHEMBL1999443.smi" -F "smarts=[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1-[#6](-[#8])=[#8])-[#8]-[#6](-[#6])=[#8]" -F "force=true" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted_forced.svg
    curl -X POST -F "file=@CHEMBL1999443.smi" -F "smarts=@aspirin.sma" -F "force=true" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.smi" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg    
    """

    number_of_files = len(request.files)
    data = None
    if number_of_files:
        if number_of_files == 1:
            data = list(request.files.values())[0].file.read()
        elif number_of_files == 2:
            data = request.files['file'].file.read()
            smarts = request.files['smarts'].file.read()
            request.params['smarts'] = smarts
    else:
        data = request.body.read()
    return highlightSmilesFragmentSvgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def highlightCtabFragmentSvgView(data, params):

    kwargs = dict()
    smarts = params.get('smarts', '')
    kwargs['size'] = int(params.get('size', 300))
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', False))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['force'] = _parseFlag(params.get('force', True))

    response.content_type = 'image/svg+xml'
    return _highlightCtabFragmentSVG(data, smarts, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/highlightCtabFragmentSvg', method=['OPTIONS', 'POST'], name="highlightCtabFragmentSvg")
def highlightCtabFragmentSvg():
    """
Converts SMILES to SVG vector graphic with a highlighted fragment described as SMARTS.
SMARTS describes the fragment to be highlighted. 
CTAB is either single molfile or SDF file. Size is the optional size of
image in pixels (default value is 300 px). 
cURL examples:

    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "computeCoords=0" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@CHEMBL1999443.mol" -F "smarts=[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1-[#6](-[#8])=[#8])-[#8]-[#6](-[#6])=[#8]" -F "force=true" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted_forced.svg
    curl -X POST -F "file=@CHEMBL1999443.mol" -F "smarts=@aspirin.sma" -F "force=true" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted_forced.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" -F "computeCoords=0" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_coords.sdf" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    """

    number_of_files = len(request.files)
    data = None
    if number_of_files:
        if number_of_files == 1:
            data = list(request.files.values())[0].file.read()
        elif number_of_files == 2:
            data = request.files['file'].file.read()
            smarts = request.files['smarts'].file.read()
            request.params['smarts'] = smarts
    else:
        data = request.body.read()
    return highlightCtabFragmentSvgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

def inchi2svgView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 300))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))

    response.content_type = 'image/svg+xml'
    return _inchi2svg(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2svg', method=['OPTIONS', 'POST'], name="inchi2svg")
def inchi2svg():
    """
Converts InChI to SVG vector graphic. This method accepts single or multiple InChIs. Size is the
optional size of image in pixels (default value is 300 px). 
cURL examples:

    curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.inchi" -F "size=400" ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
    """

    inchis = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return inchi2svgView(inchis, request.params)


# ----------------------------------------------------------------------------------------------------------------------
# Similarity maps
# ----------------------------------------------------------------------------------------------------------------------


def smiles2SimilarityMapSvgView(data, params):

    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['width'] = int(params.get('width', 500))
    kwargs['height'] = int(params.get('height', 500))
    kwargs['radius'] = int(params.get('radius', 2))
    kwargs['fingerprint'] = params.get('fingerprint', 'morgan')
    kwargs['format'] = 'svg'

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/svg+xml'
    ret = _smiles2SimilarityMapSVG(data, **kwargs)
    return ret

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2SimilarityMapSvg', method=['OPTIONS', 'POST'], name="smiles2SimilarityMapSvg")
def smiles2SimilarityMapSvg():
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires *exactly two SMILES*
cURL examples:

    curl -X POST --data-binary @sim.smi ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.svg
    curl -X POST -F "file=@sim.smi" -F "width=500" -F "height=500" -F "fingerprint=tt" ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.svg
    curl -X POST -F "file=@sim.smi" -F "width=500" -F "height=500" -F "fingerprint=ap" ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.svg
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smiles2SimilarityMapSvgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def sdf2SimilarityMapSvgView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    kwargs['width'] = int(params.get('width', 500))
    kwargs['height'] = int(params.get('height', 500))
    kwargs['radius'] = int(params.get('radius', 2))
    kwargs['fingerprint'] = params.get('fingerprint', 'morgan')
    kwargs['format'] = 'svg'

    response.content_type = 'image/svg+xml'
    ret = _sdf2SimilarityMapSVG(data, **kwargs)
    return ret

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sdf2SimilarityMapSvg', method=['OPTIONS', 'POST'], name="sdf2SimilarityMapSvg")
def sdf2SimilarityMapSvg():
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires SDF containing *exactly two mols*
cURL examples:

    curl -X POST --data-binary @sim.sdf ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.svg
    curl -X POST -F "file=@sim.sdf" -F "width=500" -F "height=500" -F "fingerprint=tt" ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.svg
    curl -X POST -F "file=@sim.sdf" -F "width=500" -F "height=500" -F "fingerprint=ap" ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.svg
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return sdf2SimilarityMapSvgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
