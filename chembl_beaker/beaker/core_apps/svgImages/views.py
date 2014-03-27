__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.core_apps.svgImages.impl import _ctab2svg, _smiles2svg
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2svg/<ctab>', method=['OPTIONS', 'GET'], name="ctab2svg")
@app.route('/ctab2svg/<ctab>/<size>', method=['OPTIONS', 'GET'], name="ctab2svg")
@app.route('/ctab2svg/<ctab>/<size>/<legend>', method=['OPTIONS', 'GET'], name="ctab2svg")
def ctab2svg(ctab, size=200, legend=''):
    """
Converts CTAB to SVG vector graphic. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles. Size is the optional size of image in pixels (default value is 200 px).
Legend is optional label in the bottom of image.
    """

    size = int(size)
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'image/svg+xml'
    return _ctab2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2svg', method=['OPTIONS', 'POST'], name="ctab2svg")
def ctab2svg():
    """
Converts CTAB to SVG vector graphic. CTAB is either single molfile or SDF file. Size is the optional size of
image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
    """

    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    legend=request.params.get('legend','')
    response.content_type = 'image/svg+xml'
    return _ctab2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2svg/<smiles>', method=['OPTIONS', 'GET'], name="smiles2svg")
@app.route('/smiles2svg/<smiles>/<size>', method=['OPTIONS', 'GET'], name="smiles2svg")
@app.route('/smiles2svg/<smiles>/<size>/<legend>', method=['OPTIONS', 'GET'], name="smiles2svg")
def smiles2svg(smiles, size=200, legend=''):
    """
Converts SMILES to SVG vector graphic. This method accepts urlsafe_base64 encoded string containing single or
multiple SMILES optionally containing header line, specific to *.smi format. Size is the optional size of image in
pixels (default value is 200 px). Legend is optional label in the bottom of image.
    """

    size = int(size)
    data = base64.urlsafe_b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    response.content_type = 'image/svg+xml'
    return _smiles2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2svg', method=['OPTIONS', 'POST'], name="smiles2svg")
def smiles2svg():
    """
Converts SMILES to SVG vector graphic. This method accepts single or multiple SMILES or *.smi file. Size is the
optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
    """

    data = request.body.read()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    size = int(request.forms.get('size', 200))
    legend=request.params.get('legend','')
    response.content_type = 'image/svg+xml'
    return _smiles2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------
