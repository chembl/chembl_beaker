__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.jsonImages.impl import _ctab2json, _smiles2json
from bottle import request, response
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2json/<ctab>', method=['OPTIONS', 'GET'], name="ctab2json")
@app.route('/ctab2json/<ctab>/<size>', method=['OPTIONS', 'GET'], name="ctab2json")
@app.route('/ctab2json/<ctab>/<size>/<legend>', method=['OPTIONS', 'GET'], name="ctab2json")
def ctab2json(ctab, size=200, legend=''):
    """
Converts CTAB to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. CTAB is urlsafe_base64 encoded string containing
single molfile or concatenation of multiple molfiles. Size is the optional size of image in pixels (default value
is 200 px). Legend is optional label in the bottom of image.
    """

    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'application/json'
    return _ctab2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2json', method=['OPTIONS', 'POST'], name="ctab2json")
def ctab2json():
    """
Converts CTAB to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. CTAB is either single molfile or SDF file. Size
is the optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
    """

    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    legend=request.params.get('legend','')
    response.content_type = 'application/json'
    return _ctab2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2json/<smiles>', method=['OPTIONS', 'GET'], name="smiles2json")
@app.route('/smiles2json/<smiles>/<size>', method=['OPTIONS', 'GET'], name="smiles2json")
@app.route('/smiles2json/<smiles>/<size>/<legend>', method=['OPTIONS', 'GET'], name="smiles2json")
def smiles2json(smiles, size=200, legend=''):
    """
Converts SMILES to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. This method accepts urlsafe_base64 encoded string
containing single or multiple SMILES optionally containing header line, specific to *.smi format. Size is the
optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
    """

    data = base64.urlsafe_b64decode(smiles)
    response.content_type = 'application/json'
    return _smiles2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2json', method=['OPTIONS', 'POST'], name="smiles2json")
def smiles2json():
    """
Converts SMILES to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. This method accepts single or multiple SMILES or
*.smi file. Size is the optional size of image in pixels (default valueis 200 px). Legend is optional label in
the bottom of image.
    """

    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    legend=request.params.get('legend','')
    response.content_type = 'application/json'
    return _smiles2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------
