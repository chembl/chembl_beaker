__author__ = 'mnowotka'

from bottle import request, response
import base64
from chembl_beaker.beaker.utils.io import _parseFlag
from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.rasterImages.impl import _ctab2image, _smiles2image

#-----------------------------------------------------------------------------------------------------------------------

def ctab2imageView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 200))
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))

    response.content_type = 'image/png'
    ret = _ctab2image(data, **kwargs)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2image/<ctab>', method=['OPTIONS', 'GET'], name="ctab2image")
def ctab2image(ctab):
    """
Converts CTAB to PNG image. CTAB is urlsafe_base64 encoded string containing
single molfile or concatenation of multiple molfiles. Size is the optional size of image in pixels (default value
is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_") > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?computeCoords=0 > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.png
    curl -X GET "${BEAKER_ROOT_URL}ctab2image/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
    curl -X GET "${BEAKER_ROOT_URL}ctab2image/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla&computeCoords=0" > out.png
    curl -X GET "${BEAKER_ROOT_URL}ctab2image/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
    curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat mcs_no_coords.sdf | base64 -w 0 | tr "+/" "-_")?legend=foo > out.png
    curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.png
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2imageView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2image', method=['OPTIONS', 'POST'], name="ctab2image")
def ctab2image():
    """
Converts CTAB to PNG image. CTAB is either single molfile or SDF file. Size is the optional size of image in pixels
(default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2image > aspirin.png
    curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
    curl -X POST -F "file=@aspirin.mol" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
    curl -X POST -F "file=@aspirin.mol" -F "legend=aspirin" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}ctab2image > out.png
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2image > out.png
    curl -X POST -F "file=@mcs_no_coords.sdf" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}ctab2image > out.png
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo" ${BEAKER_ROOT_URL}ctab2image > out.png
    curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}ctab2image > out.png
    curl -X POST -F "file=@aspirin.mol" -F "size=400" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2imageView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def smiles2imageView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 200))
    separator = params.get('separator', '|')
    kwargs['legend'] = params.get('legend', '').split(separator)
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['atomMapNumber'] = _parseFlag(params.get('atomMapNumber', False))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/png'
    ret = _smiles2image(data, **kwargs)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2image/<smiles>', method=['OPTIONS', 'GET'], name="smiles2image")
def smiles2image(smiles):
    """
Converts SMILES to PNG image. This method accepts urlsafe_base64 encoded string containing single or multiple
SMILES optionally containing header line, specific to *.smi format. Size is the optional size of image in pixels
(default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.png
    curl -X GET "${BEAKER_ROOT_URL}smiles2image/"$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
    curl -X GET "${BEAKER_ROOT_URL}smiles2image/"$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out.png
    curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out.png
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2imageView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2image', method=['OPTIONS', 'POST'], name="smiles2image")
def smiles2image():
    """
Converts SMILES to PNG image. This method accepts single or multiple SMILES or *.smi file. Size is the optional
size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST -F "file=@aspirin_no_header.smi" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST -F "file=@aspirin_with_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST -F "file=@aspirin_no_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST -F "file=@aspirin_with_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST -F "file=@aspirin_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
    curl -X POST -F "file=@mcs.smi" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}smiles2image > out.png
    curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}smiles2image > out.png
    curl -X POST -F "file=@mcs.smi" -F "legend=foo" ${BEAKER_ROOT_URL}smiles2image > out.png
    curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo" ${BEAKER_ROOT_URL}smiles2image > out.png
    curl -X POST -F "file=@mcs.smi" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > out.png
    curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > out.png
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles2imageView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------
