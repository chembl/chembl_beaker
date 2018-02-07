__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from bottle import request
import base64
import os
from chembl_beaker.beaker import app, config
from chembl_beaker.beaker.utils.io import _parseFlag
from chembl_beaker.beaker.core_apps.osra.impl import _image2ctab, _image2smiles

# ----------------------------------------------------------------------------------------------------------------------


def image2ctabView(img, params):

    kwargs = dict()
    kwargs['jaggy'] = _parseFlag(params.get('jaggy', False))
    kwargs['adaptive'] = _parseFlag(params.get('adaptive', False))
    kwargs['unpaper'] = int(params.get('unpaper', 0))

    known_location = '/usr/local/bin/osra'
    if not os.path.exists(known_location):
        known_location = '/usr/bin/osra'
    return _image2ctab(img, config.get('osra_binaries_location', known_location), **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/image2ctab/<image>', method=['OPTIONS', 'GET'], name="image2ctab")
def image2ctab(image):
    """
Uses OSRA to convert image to CTAB. Image should be urlsafe_base64 encoded data of 300 DPI png graphic.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}image2ctab/$(cat mol.jpg | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}image2ctab/$(cat mol.png | base64 -w 0 | tr "+/" "-_")
    """
    if image.startswith('data:'):
        img = base64.urlsafe_b64decode(image[image.find(',')+1:])
    else:
        img = base64.urlsafe_b64decode(image)
    return image2ctabView(img, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/image2ctab', method=['OPTIONS', 'POST'], name="image2ctab")
def image2ctab():
    """
Uses OSRA to convert image to CTAB. Image should be 300 DPI png graphic.
cURL examples:

    curl -X POST --data-binary @mol.jpg ${BEAKER_ROOT_URL}image2ctab
    curl -X POST -F "file=@mol.jpg" ${BEAKER_ROOT_URL}image2ctab
    curl -X POST --data-binary @mol.png ${BEAKER_ROOT_URL}image2ctab
    curl -X POST -F "file=@mol.png" ${BEAKER_ROOT_URL}image2ctab
    """
    img = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    if img.startswith('data:'):
        img = base64.b64decode(img[img.find(',')+1:])
    return image2ctabView(img, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def image2smilesView(img, params):

    kwargs = dict()
    kwargs['jaggy'] = _parseFlag(params.get('jaggy', False))
    kwargs['adaptive'] = _parseFlag(params.get('adaptive', False))
    kwargs['unpaper'] = int(params.get('unpaper', 0))

    known_location = '/usr/local/bin/osra'
    if not os.path.exists(known_location):
        known_location = '/usr/bin/osra'
    return _image2smiles(img, config.get('osra_binaries_location', known_location), **kwargs)


# ----------------------------------------------------------------------------------------------------------------------

@app.route('/image2smiles/<image>', method=['OPTIONS', 'GET'], name="image2smiles")
def image2smiles(image):
    """
Uses OSRA to convert image to CTAB. Image should be urlsafe_base64 encoded data of 300 DPI png graphic.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}image2smiles/$(cat mol.jpg | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}image2smiles/$(cat mol.png | base64 -w 0 | tr "+/" "-_")
    """
    if image.startswith('data:'):
        img = base64.urlsafe_b64decode(image[image.find(',')+1:])
    else:
        img = base64.urlsafe_b64decode(image)
    return image2smilesView(img, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/image2smiles', method=['OPTIONS', 'POST'], name="image2smiles")
def image2smiles():
    """
Uses OSRA to convert image to CTAB. Image should be 300 DPI png graphic.
cURL examples:

    curl -X POST --data-binary @mol.jpg ${BEAKER_ROOT_URL}image2smiles
    curl -X POST -F "file=@mol.jpg" ${BEAKER_ROOT_URL}image2smiles
    curl -X POST --data-binary @mol.png ${BEAKER_ROOT_URL}image2smiles
    curl -X POST -F "file=@mol.png" ${BEAKER_ROOT_URL}image2smiles
    """

    img = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    if img.startswith('data:'):
        img = base64.b64decode(img[img.find(',')+1:])
    return image2smilesView(img, request.params)

# ----------------------------------------------------------------------------------------------------------------------
