__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from bottle import request
import os
from beaker import app, config
from beaker.utils.io import _parseFlag
from beaker.core_apps.osra.impl import _image2ctab, _image2smiles

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
    img = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    if img.startswith(b'data:'):
        img = base64.b64decode(img[img.find(b',')+1:])
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

    img = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    if img.startswith(b'data:'):
        img = base64.b64decode(img[img.find(b',')+1:])
    return image2smilesView(img, request.params)

# ----------------------------------------------------------------------------------------------------------------------
