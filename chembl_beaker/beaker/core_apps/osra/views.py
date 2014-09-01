__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from bottle import request
import base64
import os
from chembl_beaker.beaker import app, config
from chembl_beaker.beaker.core_apps.osra.impl import _image2ctab, _image2smiles

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/image2ctab/<image>', method=['OPTIONS', 'GET'], name="image2ctab")
def image2ctab(image):
    """
Uses OSRA to convert image to CTAB. Image should be urlsafe_base65 encoded data of 300 DPI png graphic.
    """

    if image.startswith('data:'):
        img = base64.urlsafe_b64decode(image[image.find(',')+1:])
    else:
        img = base64.urlsafe_b64decode(image)
    #TODO: check /usr/local/bin/osra when no explicit value given
    known_location = '/usr/local/bin/osra'
    if not os.path.exists(known_location):
        known_location = '/usr/bin/osra'
    return _image2ctab(img, config.get('osra_binaries_location', known_location))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/image2ctab', method=['OPTIONS', 'POST'], name="image2ctab")
def image2ctab():
    """
Uses OSRA to convert image to CTAB. Image should be 300 DPI png graphic.
    """

    img = request.body.read()
    if img.startswith('data:'):
        img = base64.b64decode(img[img.find(',')+1:])
    #TODO: check /usr/local/bin/osra when no explicit value given
    known_location = '/usr/local/bin/osra'
    if not os.path.exists(known_location):
        known_location = '/usr/bin/osra'
    return _image2ctab(img, config.get('osra_binaries_location', known_location))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/image2smiles/<image>', method=['OPTIONS', 'GET'], name="image2smiles")
def image2smiles(image):
    """
Uses OSRA to convert image to CTAB. Image should be urlsafe_base65 encoded data of 300 DPI png graphic.
    """

    if image.startswith('data:'):
        img = base64.urlsafe_b64decode(image[image.find(',')+1:])
    else:
        img = base64.urlsafe_b64decode(image)
    #TODO: check /usr/local/bin/osra when no explicit value given
    known_location = '/usr/local/bin/osra'
    if not os.path.exists(known_location):
        known_location = '/usr/bin/osra'
    return _image2smiles(img, config.get('osra_binaries_location', known_location))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/image2smiles', method=['OPTIONS', 'POST'], name="image2smiles")
def image2smiles():
    """
Uses OSRA to convert image to CTAB. Image should be 300 DPI png graphic.
    """

    img = request.body.read()
    if img.startswith('data:'):
        img = base64.b64decode(img[img.find(',')+1:])
    #TODO: check /usr/local/bin/osra when no explicit value given
    known_location = '/usr/local/bin/osra'
    if not os.path.exists(known_location):
        known_location = '/usr/bin/osra'
    return _image2smiles(img, config.get('osra_binaries_location', known_location))

#-----------------------------------------------------------------------------------------------------------------------