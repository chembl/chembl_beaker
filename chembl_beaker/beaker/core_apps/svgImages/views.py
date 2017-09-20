__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.core_apps.svgImages.impl import _ctab2svg, _smiles2svg, _inchi2svg
from chembl_beaker.beaker.core_apps.svgImages.impl import _highlightCtabFragmentSVG, _highlightSmilesFragmentSVG
from chembl_beaker.beaker.utils.io import _parseFlag
import base64

# ----------------------------------------------------------------------------------------------------------------------


def ctab2svgView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 200))
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['wedgeBonds'] = _parseFlag(params.get('wedgeBonds', True))
    kwargs['fitImage'] = _parseFlag(params.get('fitImage', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))

    response.content_type = 'image/svg+xml'
    return _ctab2svg(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2svg/<ctab>', method=['OPTIONS', 'GET'], name="ctab2svg")
def ctab2svg(ctab):
    """
Converts CTAB to SVG vector graphic. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles. Size is the optional size of image in pixels (default value is 200 px).
Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_") > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?computeCoords=0 > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
    curl -X GET "${BEAKER_ROOT_URL}ctab2svg/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.svg
    curl -X GET "${BEAKER_ROOT_URL}ctab2svg/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla&computeCoords=0" > out.svg
    curl -X GET "${BEAKER_ROOT_URL}ctab2svg/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.svg
    curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat mcs_no_coords.sdf | base64 -w 0 | tr "+/" "-_")?legend=foo > out.svg    
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2svgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2svg', method=['OPTIONS', 'POST'], name="ctab2svg")
def ctab2svg():
    """
Converts CTAB to SVG vector graphic. CTAB is either single molfile or SDF file. Size is the optional size of
image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "legend=aspirin" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.mol" -F "size=400" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs_no_coords.sdf" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo" ${BEAKER_ROOT_URL}ctab2svg > out.svg
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}ctab2svg > out.svg    
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2svgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2svgView(data, params):

    kwargs = dict()
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['size'] = int(params.get('size', 200))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['wedgeBonds'] = _parseFlag(params.get('wedgeBonds', True))
    kwargs['fitImage'] = _parseFlag(params.get('fitImage', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/svg+xml'
    return _smiles2svg(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2svg/<smiles>', method=['OPTIONS', 'GET'], name="smiles2svg")
def smiles2svg(smiles):
    """
Converts SMILES to SVG vector graphic. This method accepts urlsafe_base64 encoded string containing single or
multiple SMILES optionally containing header line, specific to *.smi format. Size is the optional size of image in
pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
    curl -X GET "${BEAKER_ROOT_URL}smiles2svg/"$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.svg
    curl -X GET "${BEAKER_ROOT_URL}smiles2svg/"$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out.svg
    curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out.svg    
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2svgView(data, params=request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2svg', method=['OPTIONS', 'POST'], name="smiles2svg")
def smiles2svg():
    """
Converts SMILES to SVG vector graphic. This method accepts single or multiple SMILES or *.smi file. Size is the
optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
    curl -X POST -F "file=@mcs.smi" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs.smi" -F "legend=foo" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs.smi" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > out.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > out.svg    
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles2svgView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def highlightSmilesFragmentSvgView(data, params):

    kwargs = dict()
    smarts = params.get('smarts', '')
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['size'] = int(params.get('size', 200))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['wedgeBonds'] = _parseFlag(params.get('wedgeBonds', True))
    kwargs['fitImage'] = _parseFlag(params.get('fitImage', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['force'] = _parseFlag(params.get('force', True))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/svg+xml'
    return _highlightSmilesFragmentSVG(data, smarts, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

@app.route('/highlightSmilesFragmentSvg/<smarts>/<smiles>', method=['OPTIONS', 'GET'], name="highlightSmilesFragmentSvg")
def highlightSmilesFragmentSvg(smarts, smiles):
    """
Converts SMILES to SVG vector graphic with a highlighted fragment described as SMARTS. 
This method accepts urlsafe_base64 encoded string containing SMARTS with a fragment to be highlighted, urlsafe_base64 
encoded string containing single or multiple SMILES optionally containing header line, specific to *.smi format. 
Size is the optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(cat aspirin.sma | base64 -w 0 | tr "+/" "-_")/$(cat CHEMBL1999443.smi | base64 -w 0 | tr "+/" "-_")?force=true > aspirin_highlighted_forced.svg
    curl -X GET "${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/"$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out_highlighted.svg
    curl -X GET "${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/"$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg/$(cat aspirin.sma | base64 -w 0 | tr "+/" "-_")/$(cat CHEMBL1999443.smi | base64 -w 0 | tr "+/" "-_")?force=true > out_highlighted_forced.svg    
    """

    data = base64.urlsafe_b64decode(smiles)

    params = request.params
    params['smarts'] = base64.urlsafe_b64decode(smarts)

    return highlightSmilesFragmentSvgView(data, params=params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/highlightSmilesFragmentSvg', method=['OPTIONS', 'POST'], name="highlightSmilesFragmentSvg")
def highlightSmilesFragmentSvg():
    """
Converts SMILES to SVG vector graphic with a highlighted fragment described as SMARTS. 
This method accepts SMARTS and single or multiple SMILES or *.smi file. Size is the
optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" -F "legend=aspirin" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_no_header.smi" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "smarts=c1ccccc1" -F "legend=aspirin" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin_with_header.smi" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@CHEMBL1999443.smi" -F "smarts=[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1-[#6](-[#8])=[#8])-[#8]-[#6](-[#6])=[#8]" -F "force=true" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted_forced.svg
    curl -X POST -F "file=@CHEMBL1999443.smi" -F "smarts=@aspirin.sma" -F "force=true" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.smi" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.smi" -F "smarts=c1ccccc1" -F "legend=foo" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "smarts=c1ccccc1" -F "legend=foo" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.smi" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_header.smi" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}highlightSmilesFragmentSvg > out_highlighted.svg    
    """

    number_of_files = len(request.files)
    data = None
    if number_of_files:
        if number_of_files == 1:
            data = request.files.values()[0].file.read()
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
    kwargs['size'] = int(params.get('size', 200))
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['wedgeBonds'] = _parseFlag(params.get('wedgeBonds', True))
    kwargs['fitImage'] = _parseFlag(params.get('fitImage', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['force'] = _parseFlag(params.get('force', True))

    response.content_type = 'image/svg+xml'
    return _highlightCtabFragmentSVG(data, smarts, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/highlightCtabFragmentSvg/<smarts>/<ctab>', method=['OPTIONS', 'GET'], name="highlightCtabFragmentSvg")
def highlightCtabFragmentSvg(smarts, ctab):
    """
Converts SMILES to SVG vector graphic with a highlighted fragment described as SMARTS. 
SMARTS is urlsafe_base64 encoded string containing fragment ot ge highlighted.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles. 
Size is the optional size of image in pixels (default value is 200 px).
Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_") > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?computeCoords=0 > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(cat aspirin.sma | base64 -w 0 | tr "+/" "-_")/$(cat CHEMBL1999443.mol | base64 -w 0 | tr "+/" "-_")?force=true > out_highlighted_forced.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out_highlighted.svg
    curl -X GET ${BEAKER_ROOT_URL}highlightCtabFragmentSvg/$(echo c1ccccc1 | base64 -w 0 | tr "+/" "-_")/$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla&computeCoords=0" > out_highlighted.svg    
    """

    data = base64.urlsafe_b64decode(ctab)

    params = request.params
    params['smarts'] = base64.urlsafe_b64decode(smarts)

    return highlightCtabFragmentSvgView(data, params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/highlightCtabFragmentSvg', method=['OPTIONS', 'POST'], name="highlightCtabFragmentSvg")
def highlightCtabFragmentSvg():
    """
Converts SMILES to SVG vector graphic with a highlighted fragment described as SMARTS.
SMARTS describes the fragment to be highlighted. 
CTAB is either single molfile or SDF file. Size is the optional size of
image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "computeCoords=0" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "legend=aspirin" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@aspirin.mol" -F "smarts=c1ccccc1" -F "size=400" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > aspirin_highlighted.svg
    curl -X POST -F "file=@CHEMBL1999443.mol" -F "smarts=[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1-[#6](-[#8])=[#8])-[#8]-[#6](-[#6])=[#8]" -F "force=true" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted_forced.svg
    curl -X POST -F "file=@CHEMBL1999443.mol" -F "smarts=@aspirin.sma" -F "force=true" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted_forced.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" -F "computeCoords=0" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs_no_coords.sdf" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" -F "legend=foo" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    curl -X POST -F "file=@mcs.sdf" -F "smarts=c1ccccc1" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}highlightCtabFragmentSvg > out_highlighted.svg
    """

    number_of_files = len(request.files)
    data = None
    if number_of_files:
        if number_of_files == 1:
            data = request.files.values()[0].file.read()
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
    kwargs['size'] = int(params.get('size', 200))
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['kekulize'] = _parseFlag(params.get('kekulize', True))
    kwargs['wedgeBonds'] = _parseFlag(params.get('wedgeBonds', True))
    kwargs['fitImage'] = _parseFlag(params.get('fitImage', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))

    response.content_type = 'image/svg+xml'
    return _inchi2svg(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2svg/<inchi>', method=['OPTIONS', 'GET'], name="inchi2svg")
def inchi2svg(inchi):
    """
Converts InChI to SVG vector graphic. This method accepts urlsafe_base64 encoded string containing single or
multiple InChIs. Size is the optional size of image in pixels (default value is 200 px).
Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}inchi2svg/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_") > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}inchi2svg/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
    curl -X GET ${BEAKER_ROOT_URL}inchi2svg/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
    """

    inchis = base64.urlsafe_b64decode(inchi)
    return inchi2svgView(inchis, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2svg', method=['OPTIONS', 'POST'], name="inchi2svg")
def inchi2svg():
    """
Converts InChI to SVG vector graphic. This method accepts single or multiple InChIs. Size is the
optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.inchi" -F "legend=aspirin" ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
    curl -X POST -F "file=@aspirin.inchi" -F "size=400" ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
    """

    inchis = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return inchi2svgView(inchis, request.params)

# ----------------------------------------------------------------------------------------------------------------------

