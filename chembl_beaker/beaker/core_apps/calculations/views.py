__author__ = 'mnowotka'

from chembl_beaker.beaker import app
from bottle import request
from chembl_beaker.beaker.core_apps.calculations.impl import _kekulize
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/kekulize/<ctab>', method=['OPTIONS', 'GET'], name="kekulize")
def kekulize(ctab):
    """
Performs kekulisation on input compounds. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
    """

    data = base64.urlsafe_b64decode(ctab)
    return _kekulize(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/kekulize', method=['OPTIONS', 'POST'], name="kekulize")
def kekulize():
    """
Performs kekulisation on input compounds. CTAB is either single molfile or SDF file.
    """

    data = request.body.read()
    return _kekulize(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sanitize/<ctab>', name="sanitize")
def sanitize(ctab):
    """
This method is not implemented yet.
    """

    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/atomIsInRing/<ctab>/<index>/<size>', name="atomIsInRing")
def atomIsInRing(ctab, index, size):
    """
This method is not implemented yet.
    """

    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/symmSSSR/<ctab>', name="symmSSSR")
def symmSSSR(ctab):
    """
This method is not implemented yet.
    """

    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sssr/<ctab>', name="sssr")
def sssr(ctab):
    """
This method is not implemented yet.
    """

    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/addHs/<ctab>', name="addHs")
def addHs(ctab):
    """
This method is not implemented yet.
    """

    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/removeHs/<ctab>', name="removeHs")
def removeHs(ctab):
    """
This method is not implemented yet.
    """

    pass

#-----------------------------------------------------------------------------------------------------------------------
