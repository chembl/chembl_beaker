__author__ = 'mnowotka'

from bottle import request, response
import base64
from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.rasterImages.impl import _ctab2image, _smiles2image

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2image/<ctab>', method=['OPTIONS', 'GET'], name="ctab2image")
@app.route('/ctab2image/<ctab>/<size>', method=['OPTIONS', 'GET'], name="ctab2image")
@app.route('/ctab2image/<ctab>/<size>/<legend>', method=['OPTIONS', 'GET'], name="ctab2image")
def ctab2image(ctab, size=200, legend=''):
    """
Converts CTAB to PNG image. CTAB is urlsafe_base64 encoded string containing
single molfile or concatenation of multiple molfiles. Size is the optional size of image in pixels (default value
is 200 px). Legend is optional label in the bottom of image.
    """

    size = int(size)
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'image/png'
    ret = _ctab2image(data,size,legend)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2image', method=['OPTIONS', 'POST'], name="ctab2image")
def ctab2image():
    """
Converts CTAB to PNG image. CTAB is either single molfile or SDF file. Size is the optional size of image in pixels
(default value is 200 px). Legend is optional label in the bottom of image.
    """

    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    legend=request.params.get('legend','')
    response.content_type = 'image/png'
    ret = _ctab2image(data,size,legend)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2image/<smiles>', method=['OPTIONS', 'GET'], name="smiles2image")
@app.route('/smiles2image/<smiles>/<size>', method=['OPTIONS', 'GET'], name="smiles2image")
@app.route('/smiles2image/<smiles>/<size>/<legend>', method=['OPTIONS', 'GET'], name="smiles2image")
def smiles2image(smiles, size=200, legend=''):
    """
Converts SMILES to PNG image. This method accepts urlsafe_base64 encoded string containing single or multiple
SMILES optionally containing header line, specific to *.smi format. Size is the optional size of image in pixels
(default value is 200 px). Legend is optional label in the bottom of image.
    """

    size = int(size)
    data = base64.urlsafe_b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    response.content_type = 'image/png'
    ret = _smiles2image(data,size,legend)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2image', method=['OPTIONS', 'POST'], name="smiles2image")
def smiles2image():
    """
Converts SMILES to PNG image. This method accepts single or multiple SMILES or *.smi file. Size is the optional
size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
    """

    data = request.body.read()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    size = int(request.forms.get('size', 200))
    legend=request.params.get('legend','')
    response.content_type = 'image/png'
    ret = _smiles2image(data,size,legend)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------
